import os
import networkx as nx
import cv2
from PIL import Image
import numpy as np
from skimage.morphology import skeletonize
from src.function_trace_source import trace_source
from src.function_label_correct import pass_label
from src.function_topolpgy_adjustment import topology_adjustment1, topology_adjustment2, topology_adjustment3
from src.mul_ske import mul_ske_tograph
from src.MulSke_graph_optimize import optimize_graph
from src.shorted_path_dijkstra import find_shorted_path
from src.function_pruning import Pruning, Connecting1, Connecting0
from src.Nearest_neighbor_label_passing import nearest_neighbor_label_passing


def main(args):
    
    image = Image.open(args.pred_v_path)
    img = cv2.imread(args.fundus_img_path)
    image_label0 = Image.open(args.pred_av_path)
    imagelabel_array = np.array(image_label0)
    h, w = imagelabel_array.shape
    image_label = np.zeros(shape=(h, w, 3))
    for i in range(0, h):
        for j in range(0, w):
            if imagelabel_array[i, j] == 2:
                image_label[i, j, :] = [255, 0, 0]
            if imagelabel_array[i, j] == 3:
                image_label[i, j, :] = [0, 0, 255]
            if imagelabel_array[i, j] == 1:
                image_label[i, j, :] = [0, 255, 0]
    image_zeros = np.zeros([h, w])
    image = np.array(image)
    img_lab = cv2.cvtColor(img, cv2.COLOR_BGR2LAB)

    # -----------------------------------------to get Vascular topology----------------------------------------------------
    i = 0
    for i in range(255, 80, -25):
        ret, st = cv2.threshold(image, i, 255, cv2.THRESH_BINARY)
        st = st / 255
        img_sk = skeletonize(st)
        image_zeros[np.array(img_sk) != 0] = 255
    ret, st = cv2.threshold(image, i, 255, cv2.THRESH_BINARY)
    st = st / 255
    img_sk1 = skeletonize(st)
    img_sk1 = img_sk1 * 255
    img_gray = st
    image_zeros_copy = image_zeros.copy()

    g, B_Gu = mul_ske_tograph(image_zeros=image_zeros, h=h, w=w)

    G_uc1, D_t, image_zeros_205, node_remove, B_intensive_G_uc, G_uc_copy, G_uc1_copy = optimize_graph(h=h, w=w, img_sk1=img_sk1, g=g, B_Gu=B_Gu, image_zeros=image_zeros, image=image,
                        image_zeros_copy=image_zeros_copy, image_label=image_label)

    G_uc = G_uc1

    b_G_uc = []
    b1_G_uc = []
    for node in G_uc:
        if G_uc.degree[node] > 2:
            b_G_uc.append(node)
        if G_uc.degree[node] == 1:
            b1_G_uc.append(node)

    G_uc = topology_adjustment3(G_uc)
    
    D_t = []
    for node in G_uc:
        if G_uc.degree[node] == 1:
            D_t.append(node)

    D_t_1 = []
    for node in G_uc:
        if G_uc.degree[node] == 1:
            D_t_1.append(node)

    # D_new = []
    for node in D_t_1:
        if node not in b1_G_uc:
            # D_new.append(node)
            D_t.append(node)


    D_G_uc1_copy = []
    for node in G_uc1_copy:
        if G_uc1_copy.degree[node] == 1:
            D_G_uc1_copy.append(node)

    D_t_x = []
    for node in b1_G_uc:
        if node not in D_G_uc1_copy:
            D_t_x.append(node)

    g0 = nx.Graph()
    H = args.H  # the optic disc center (H).
    g3, gl = find_shorted_path(G_uc=G_uc, img_lab=img_lab, img_gray=img_gray, image=image, H=H, D_t=D_t, g0=g0)
    bn_g3 = []
    for node in g3:
        if g3.degree[node] > 2:
            bn_g3.append(node)

    b_gl = []
    D_gl = []
    for node in gl:
        if gl.degree[node] > 2:
            b_gl.append(node)
        if gl.degree[node] == 1:
            D_gl.append(node)

    for c in sorted(nx.connected_components(g3), key=len, reverse=True):  
        c = g3.subgraph(c).copy()
        if len(c.nodes) <= 4:
            for node in c.nodes:
                g3.remove_node(node)

    g4 = topology_adjustment1(g3)

    b_g4_1 = []
    for node in g4.nodes:
        if g4.degree[node] > 2:
            b_g4_1.append(node)

    g4 = topology_adjustment2(g4)

    b4_g4 = []
    b3_g4 = []
    b1_g4 = []
    for node in g4.nodes:
        if g4.degree[node] == 4:
            b4_g4.append(node)
        if g4.degree[node] == 3:
            b3_g4.append(node)
        if g4.degree[node] == 1:
            b1_g4.append(node)

    image_circle = np.zeros(shape=image.shape)
    K_s = args.Ks
    Disc_size = int(25 * K_s)
    cv2.circle(image_circle, (H[1], H[0]), Disc_size, (128, 0, 0), -1)
    Disk_node = []
    for i in range(0, image_circle.shape[0]):
        for j in range(0, image_circle.shape[1]):
            if image_circle[i, j] == 128:
                Disk_node.append((i, j))

    g4 = Pruning(g4=g4, Disk_node=Disk_node, image_label=image_label, ratio=K_s)
    b4_g4 = []
    b3_g4 = []
    b1_g4 = []
    for node in g4.nodes:
        if g4.degree[node] == 4:
            b4_g4.append(node)
        if g4.degree[node] == 3:
            b3_g4.append(node)
        if g4.degree[node] == 1:
            b1_g4.append(node)

    g4 = Connecting0(g4=g4, Disk_node=Disk_node, Disc_size=Disc_size, H=H, K=2, K_s=K_s)
    b4_g4 = []
    b3_g4 = []
    b1_g4 = []
    for node in g4.nodes:
        if g4.degree[node] == 4:
            b4_g4.append(node)
        if g4.degree[node] == 3:
            b3_g4.append(node)
        if g4.degree[node] == 1:
            b1_g4.append(node)

    # -------------------------------------- Vascular Directed Topology Estimation------------------------------------------
    G_direct, G_direct_old, g4_not_cut, g4, close_disc = trace_source(g4=g4, Disk_node=Disk_node, H=H, Disc_size=Disc_size, D_t_x = D_t_x, image_label=image_label)

    G_direct_new = G_direct.reverse()

    starting_point = []
    for node in G_direct_new:
        in_degree = (
            0 if len({G_direct_new.in_degree(nbunch=node)}) == 0 else G_direct_new.in_degree(nbunch=node))
        out_degree = (
            0 if len({G_direct_new.out_degree(nbunch=node)}) == 0 else G_direct_new.out_degree(nbunch=node))
        if in_degree == 0 and out_degree != 0:
            starting_point.append(node)

    G_direct_new, g4, g4_not_cut, starting_point = Connecting1(G_direct_new, g4, g4_not_cut, starting_point, Disc_size, H, close_disc, K=1)

    # --------------------------------------------------------A/V label delivery--------------------------------------------
    G_direct_v, G_direct_a, G_direct_twoclasses, G_direct_v_allseg, G_direct_a_allseg, G_direct_allseg = pass_label(starting_point, G_direct_new, g4, g4_not_cut, image_label, h, H, Disc_size)
    BN_5_ = []
    BN_4_ = []
    BN_3_ = []
    for node in G_direct_twoclasses.nodes:
        if G_direct_twoclasses.degree[node] > 4:
            BN_5_.append(node)
        if G_direct_twoclasses.degree[node] == 4:
            BN_4_.append(node)
        if G_direct_twoclasses.degree[node] == 3:
            BN_3_.append(node)

    nrb_out = []
    for edge in list(G_direct_twoclasses.out_edges(nbunch=(237, 184))):
        nrb_out.append(edge[1])
    node_test = []
    node_test0 = []
    node_test1 = []
    for node in BN_3_:
        in_degree = (
            0 if len({G_direct_v.in_degree(nbunch=node)}) == 0 else G_direct_v.in_degree(nbunch=node))
        out_degree = (
            0 if len({G_direct_v.out_degree(nbunch=node)}) == 0 else G_direct_v.out_degree(nbunch=node))
        if in_degree == 1 and out_degree == 2 and node not in G_direct_a.nodes:
            node_test.append(node)
        if in_degree == 1 and out_degree == 2 and node in G_direct_v.nodes:
            node_test0.append(node)
        if in_degree == 1 and out_degree == 1:
            node_test1.append(node)

    # ---------------------------------------optimization of arteriovenous classification-----------------------------------
    G_direct_v_allseg_copy = G_direct_v_allseg.copy()
    G_direct_a_allseg_copy = G_direct_a_allseg.copy()
    for node in G_direct_v_allseg_copy:
        if G_direct_v_allseg.degree[node] == 0:
            G_direct_v_allseg.remove_node(node)
            G_direct_allseg.remove_node(node)
    for node in G_direct_a_allseg_copy:
        if G_direct_a_allseg.degree[node] == 0:
            G_direct_a_allseg.remove_node(node)
            G_direct_allseg.remove_node(node)

    image = nearest_neighbor_label_passing(h=h, w=w, image=image_label, G_direct_twoclasses=G_direct_allseg,
                                        G_direct_v=G_direct_v_allseg, G_direct_a=G_direct_a_allseg, Disk_node=Disk_node,
                                        g4=g4_not_cut)

    save_path = os.path.join(args.opti_save_path, 'image_vessel_correct.png')
    cv2.imwrite(save_path, image)

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Artery-Vein-Optimization")

    # parser.add_argument("--fundus_img_path", help="the path of color fundus image")
    # parser.add_argument("--pred_av_path", type=str, help="the path of perdition artery and vein result")
    # parser.add_argument("--pred_v_path", type=str, help="the path of perdition vessel result")
    # parser.add_argument("--opti_save_path", type=str, help="the save path of optimization result")
    # parser.add_argument("--H", help="Center point coordinates of the optic disc")
    # parser.add_argument("--Ks", type=float, help=" the value of Ks determined by the size of the image, hyperparameter")

    # parser.add_argument("--fundus_img_path", default="./data/AV-DRIVE/orig/1_image.png", help="the path of color fundus image")
    # parser.add_argument("--pred_av_path", default="./data/AV-DRIVE/pred_result_cnn/1_pred_label.png", type=str, help="the path of perdition artery and vein result")
    # parser.add_argument("--pred_v_path", default="./data/AV-DRIVE/pred_result_cnn/1_pred_vessel.png", type=str, help="the path of perdition vessel result")
    # parser.add_argument("--opti_save_path", default="./data/AV-DRIVE", type=str, help="the save path of optimization result")
    # parser.add_argument("--H", default=(290, 485), help="Center point coordinates of the optic disc")
    # parser.add_argument("--Ks", default=1.0, type=float, help=" the value of Ks determined by the size of the image, hyperparameter")

    # parser.add_argument("--fundus_img_path", default="./data/LES/orig/3_image.png", help="the path of color fundus image")
    # parser.add_argument("--pred_av_path", default="./data/LES/pred_result_cnn/3_pred_label.png", type=str, help="the path of perdition artery and vein result")
    # parser.add_argument("--pred_v_path", default="./data/LES/pred_result_cnn/3_pred_vessel.png", type=str, help="the path of perdition vessel result")
    # parser.add_argument("--opti_save_path", default="./data/LES", type=str, help="the save path of optimization result")
    # parser.add_argument("--H", default=(720, 803), help="Center point coordinates of the optic disc")
    # parser.add_argument("--Ks", default=2.5, type=float, help=" the value of Ks determined by the size of the image, hyperparameter")

    parser.add_argument("--fundus_img_path", default="./data/CY-BIT/orig/5126_image.png", help="the path of color fundus image")
    parser.add_argument("--pred_av_path", default="./data/CY-BIT/pred_result_cnn/5126_pred_label.png", type=str, help="the path of perdition artery and vein result")
    parser.add_argument("--pred_v_path", default="./data/CY-BIT/pred_result_cnn/5126_pred_vessel.png", type=str, help="the path of perdition vessel result")
    parser.add_argument("--opti_save_path", default="./data/CY-BIT", type=str, help="the save path of optimization result")
    parser.add_argument("--H", default=(745, 751), help="Center point coordinates of the optic disc")
    parser.add_argument("--Ks", default=2.5, type=float, help=" the value of Ks determined by the size of the image, hyperparameter")

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = parse_args()

    main(args)