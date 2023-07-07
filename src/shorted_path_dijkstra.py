import networkx as nx
import cv2
import numpy as np
from src.color_diff import delta_e_cie2000
from colormath.color_objects import LabColor
from cv2 import DIST_L2, DIST_MASK_3
import math

def edge_weight(n1, n2, img1, img2, img3, w_v=1, w_c=1, w_w=1):
    color1 = img1[n1]
    color2 = img1[n2]
    color1 = LabColor(lab_l=color1[0], lab_a=color1[1], lab_b=color1[2])
    color2 = LabColor(lab_l=color2[0], lab_a=color2[1], lab_b=color2[2])
    c_ = delta_e_cie2000(color1=color1, color2=color2, Kl=1, Kc=1, Kh=1)
    bw_dist = cv2.distanceTransform(src=img2.astype(np.uint8), distanceType=DIST_L2, maskSize=DIST_MASK_3)
    w_ = bw_dist[n1] + bw_dist[n2]
    p_ = 255 - ((int(img3[n1]) + int(img3[n2])) / 2)
    w_n12n2 = math.exp((w_v * math.log(1 + p_)) + (w_c * math.log(1 + c_)) + (w_w * math.log(1 + w_)))
    return w_n12n2

def edge_weight2(n1, n2, img1, img2, img3, w1=1.0, w2=1.0, w3=1.0):
    color1 = img1[n1]
    color2 = img1[n2]
    color1 = LabColor(lab_l=color1[0], lab_a=color1[1], lab_b=color1[2])
    color2 = LabColor(lab_l=color2[0], lab_a=color2[1], lab_b=color2[2])
    c_ = delta_e_cie2000(color1=color1, color2=color2, Kl=1, Kc=1, Kh=1)
    bw_dist = cv2.distanceTransform(src=img2.astype(np.uint8), distanceType=DIST_L2, maskSize=DIST_MASK_3)
    w_ = bw_dist[n1] + bw_dist[n2]
    p_ = 255 - ((int(img3[n1]) + int(img3[n2])) / 2)
    w_n12n2 = (w3 * p_) / (math.exp((w1 * w_ + w2 * c_)))
    return w_n12n2

def edge_weight1(n1, n2, img2, img3, w=1.0):
    bw_dist = cv2.distanceTransform(src=img2.astype(np.uint8), distanceType=DIST_L2,
                                    maskSize=DIST_MASK_3)
    w_ = bw_dist[n1] + bw_dist[n2]
    p_ = 255 - ((int(img3[n1]) + int(img3[n2])) / 2)
    w_n12n2 = p_ / (math.exp((w * w_)))
    return w_n12n2


def find_shorted_path(G_uc, img_lab, img_gray, image, H, D_t, g0):
    for edge in G_uc.edges:
        (n1, n2) = edge
        w_n12n2 = edge_weight(n1, n2, img_lab, img_gray, image, w_v=2, w_c=2, w_w=1)
        G_uc.add_weighted_edges_from([(n1, n2, w_n12n2)])

    m = 0
    Dt_cg1 = []
    for c in sorted(nx.connected_components(G_uc), key=len, reverse=True):
        D_c = []
        Dt_cg = []
        c = G_uc.subgraph(c).copy()
        if len(c.nodes) == 1:
            break
        c_g = nx.Graph()
        src = H
        for node in c:
            c_g.add_node(node)
            degree = G_uc.degree[node]
            if degree == 1:
                D_c.append(node)
        if len(D_c) == 0:
            break
        if H not in c:
            src = D_c[0]
        for node in D_t:
            if node in c_g:
                Dt_cg.append(node)
                Dt_cg1.append(node)
        for t1 in Dt_cg:
            if src != t1:
                path = nx.dijkstra_path(c, source=src,
                                        target=t1)
                nx.add_path(g0, path)

    g_ = nx.Graph()
    g1 = g0.copy()
    gw = G_uc.copy()
    gl = gw.copy()
    for node in g1.nodes:
        gl.remove_node(node)

    for c in sorted(nx.connected_components(gl), key=len, reverse=True):
        D_gl = []
        c = gl.subgraph(c).copy()
        if len(c.nodes) < 4:
            continue
        for edge in c.edges:
            (n1, n2) = edge
            w_n12n2 = (edge_weight1(n1, n2, img_gray, image, w=5))
            G_uc.add_weighted_edges_from([(n1, n2, w_n12n2)])

    g2 = g1.copy()
    for c in sorted(nx.connected_components(G_uc), key=len, reverse=True):
        D_c = []
        Dt_cg = []
        c = G_uc.subgraph(c).copy()
        if len(c.nodes) == 1:
            break
        c_g = nx.Graph()
        src = H
        for node in c:
            c_g.add_node(node)
            degree = G_uc.degree[node]
            if degree == 1:
                D_c.append(node)
        if len(D_c) == 0:
            break
        if H not in c:
            src = D_c[0]
        for node in D_t:
            if node in c_g:
                Dt_cg.append(node)
        for t1 in Dt_cg:
            if src != t1:
                path = nx.dijkstra_path(c, source=src, target=t1)
                nx.add_path(g2, path)
                nx.add_path(g_, path)

    g3 = g2.copy()
    gl = G_uc.copy()
    for node in g3.nodes:
        gl.remove_node(node)
    for c in sorted(nx.connected_components(gl), key=len, reverse=True):
        D_gl = []
        c = gl.subgraph(c).copy()
        normal_node = []
        for node in c.nodes:
            if G_uc.degree[node] == 2:
                normal_node.append(node)
        if len(normal_node) < 2:
            continue
        m = 0
        for node in c.nodes:
            if G_uc.degree[node] > 2:
                m += 1
        if m / len(c.nodes) > 0.4:
            continue
        for edge in c.edges:
            (n1, n2) = edge
            w_n12n2 = (edge_weight1(n1, n2, img_gray, image, w=5))
            G_uc.add_weighted_edges_from([(n1, n2, w_n12n2)])
        for node in c.nodes:
            if c.degree[node] == 1:
                D_gl.append(node)
        if len(D_gl) == 0:
            if len(c.nodes) < 10:
                continue
            else:
                node_near_g3 = []
                for node in c:
                    for edge in G_uc.edges(nbunch=node):
                        if edge[1] in g3 and node not in node_near_g3:
                            node_near_g3.append(node)
                if len(node_near_g3) != 0:
                    D_gl.append(node_near_g3[-1])
                if len(D_gl) == 1:
                    node_near_g3 = []
                    bn_3 = []
                    for node in c:
                        if G_uc.degree[node] == 3:
                            bn_3.append(node)
                    for node in bn_3:
                        nrb = []
                        for edge in list(G_uc.edges(nbunch=node)):
                            nrb.append(edge[1])
                        for node1 in nrb:
                            if G_uc.degree[node1] == 2:
                                for edge1 in list(G_uc.edges(nbunch=node1)):
                                    if edge1[1] != node and edge1[1] in nrb:
                                        node_near_g3.append(node1)
                if len(node_near_g3) != 0:
                    D_gl.append(node_near_g3[-1])
        if len(D_gl) == 0:
            continue

        if len(D_gl) == 1:
            node_near_g3 = []
            for node in c:
                for edge in G_uc.edges(nbunch=node):
                    if edge[1] in g3 and node not in node_near_g3 and node != D_gl[0]:
                        node_near_g3.append(node)
            if len(node_near_g3) == 0:
                bn_3 = []
                for node in c:
                    if G_uc.degree[node] == 3:
                        bn_3.append(node)
                for node in bn_3:
                    nrb = []
                    for edge in list(G_uc.edges(nbunch=node)):
                        nrb.append(edge[1])
                    for node1 in nrb:
                        if G_uc.degree[node1] == 2:
                            for edge1 in list(G_uc.edges(nbunch=node1)):
                                if edge1[1] != node and edge1[1] in nrb:
                                    node_near_g3.append(node1)
            if len(node_near_g3) != 0:
                D_gl.append(node_near_g3[-1])

        src = D_gl.pop(-1)
        for t1 in D_gl:
            path = [(0, 0)]
            if src != t1:
                path = nx.dijkstra_path(c, source=src, target=t1)
                nx.add_path(g3, path)
            other_endpoint_path = path[0]
            if gw.degree[t1] != 1:
                for i in range(0, len(gw.edges(nbunch=t1))):
                    q = list(gw.edges(nbunch=t1))[i][1]
                    if q not in c.nodes:
                        g3.add_edge(t1, q)
            if gw.degree[other_endpoint_path] != 1:
                for i in range(0, len(gw.edges(nbunch=other_endpoint_path))):
                    q = list(gw.edges(nbunch=other_endpoint_path))[i][1]
                    if q not in c.nodes:
                        g3.add_edge(other_endpoint_path, q)

    return g3, gl