import networkx as nx
import cv2
import numpy as np
from skimage.morphology import skeletonize
from networkx.algorithms import components
def neighborhood_1(i, j, image, h, k1=0, k2=3):
    count = 0
    if i == 0:
        k1, k2 = 1, 3
    if i == h - 1:
        k1, k2 = 0, 2
    for hang in range(k1, k2):
        for lie in range(0, 3):
            if hang == lie == 1:
                continue
            if image[i - 1 + hang, j - 1 + lie] == 255:
                count += 1
    return count


def Neighborhood_8(i, j, image, h, k1=0, k2=3):
    count = 0
    neighborhood_8 = []
    if i == 0:
        k1, k2 = 1, 3
    if i == h - 1:
        k1, k2 = 0, 2
    for hang in range(k1, k2):
        for lie in range(0, 3):
            if hang == lie == 1:
                continue
            if image[i - 1 + hang, j - 1 + lie] == 255:
                count += 1
                neighborhood_8.append((i - 1 + hang, j - 1 + lie))
    return count, neighborhood_8

def find_10x10(i, j, image):
    count = 0
    B_intensive_G_uc_sub = []
    for m in range(-5, 5):
        for n in range(-5, 5):
            if image[i + m, j + n] == 1:
                B_intensive_G_uc_sub.append((i + m, j + n))
                count += 1
    return count, B_intensive_G_uc_sub


def find_16x16(i, j, image):
    count = 0
    B_intensive_G_uc_sub = []
    for m in range(-8, 8):
        for n in range(-8, 8):
            if image[i + m, j + n] == 1:
                B_intensive_G_uc_sub.append((i + m, j + n))
                count += 1
    return count, B_intensive_G_uc_sub

def find_16x16_new(i, j, image):
    count = 0
    B_intensive_G_uc_sub = []
    for m in range(-3, 13):
        for n in range(-3, 13):
            if image[i + m, j + n] == 1:
                B_intensive_G_uc_sub.append((i + m, j + n))
                count += 1
    return count, B_intensive_G_uc_sub

def optimize_1th(B_Gu, G_uc, image_zeros, D_t, h):
    while len(B_Gu) > 0:
        b = B_Gu.pop(-1)
        if b not in G_uc:
            continue
        (i, j) = b
        count, neighborhood_8_b = Neighborhood_8(i=i, j=j, image=image_zeros, h=h, k1=0, k2=3)
        for item, n in enumerate(neighborhood_8_b):
            if n in D_t:
                continue
            if n not in G_uc.nodes:
                continue
            (i, j) = n
            count, neighborhood_8_n = Neighborhood_8(i, j, image_zeros, h, 0, 3)
            for _, m in enumerate(neighborhood_8_n):
                if m == b or m in neighborhood_8_b:
                    continue
                if m not in G_uc.nodes:
                    continue
                G_uc.add_edge(b, m)
            if n in G_uc:
                (i, j) = n
                image_zeros[i, j] = 0
                count1, neighborhood_8_n = Neighborhood_8(i, j, image_zeros, h, 0, 3)
                if G_uc.degree[n] > count1:
                    for i in range(0, len(G_uc.edges(nbunch=n))):
                        q = list(G_uc.edges(nbunch=n))[i][1]
                        if abs(q[0] - n[0]) > 1 or abs(q[1] - n[1]) > 1:
                            G_uc.add_edge(b, q)
                G_uc.remove_node(n)
    return G_uc

def optimize_2th(G_uc1, image_zeros_new, D_t, h):
    B_G_uc1 = []
    for node in G_uc1.nodes:
        if G_uc1.degree[node] > 2:
            B_G_uc1.append(node)
    while len(B_G_uc1) > 0:
        b = B_G_uc1.pop(-1)
        if b not in G_uc1:
            continue
        (i, j) = b
        count, neighborhood_8_b = Neighborhood_8(i=i, j=j, image=image_zeros_new, h=h,  k1=0, k2=3)
        for item, n in enumerate(neighborhood_8_b):
            if n in D_t:
                continue
            if n not in G_uc1.nodes:
                continue
            (i, j) = n
            count, neighborhood_8_n = Neighborhood_8(i, j, image_zeros_new, h, 0, 3)
            for item, m in enumerate(neighborhood_8_n):
                if m == b or m in neighborhood_8_b:
                    continue
                if m not in G_uc1.nodes:
                    continue
                G_uc1.add_edge(b, m)
            if n in G_uc1:
                (i, j) = n
                image_zeros_new[i, j] = 0
                count1, neighborhood_8_n = Neighborhood_8(i, j, image_zeros_new, h, 0, 3)
                if G_uc1.degree[n] > count1:
                    for i in range(0, len(G_uc1.edges(nbunch=n))):
                        q = list(G_uc1.edges(nbunch=n))[i][1]
                        if abs(q[0] - n[0]) > 1 or abs(q[1] - n[1]) > 1:
                            G_uc1.add_edge(b, q)
                G_uc1.remove_node(n)
    return G_uc1

def optimize_graph(h, w, img_sk1, g, B_Gu, image_zeros, image, image_zeros_copy, image_label):
    D_t = []
    for i in range(0, h):
        for j in range(0, w):
            if img_sk1[i, j] != 0:
                g.add_node((i, j))
                k1, k2 = 0, 3
                count = neighborhood_1(i=i, j=j, image=img_sk1, h=h, k1=k1, k2=k2)
                if count == 1:
                    D_t.append((i, j))

    G_uc = g.copy()
    G_uc1 = g.copy()
    B_intensive_G_uc1 = []
    for node in G_uc1:
        if G_uc1.degree[node] > 4:
            B_intensive_G_uc1.append(node)
    img_B_intensive_G_uc1 = np.zeros([h, w])
    for node in B_intensive_G_uc1:
        img_B_intensive_G_uc1[node[0], node[1]] = 1

    G_uc = optimize_1th(B_Gu=B_Gu, G_uc=G_uc, image_zeros=image_zeros, D_t=D_t, h=h)

    B_intensive_G_uc = []
    for node in G_uc:
        if G_uc.degree[node] > 2:
            B_intensive_G_uc.append(node)
    ret1, st1 = cv2.threshold(image, 180, 255, cv2.THRESH_BINARY)
    st1 = st1 / 255
    img_sk205 = skeletonize(st1)
    img_sk205 = img_sk205 * 255
    image_zeros_205 = np.zeros([h, w])
    image_zeros_205[np.array(img_sk205) != 0] = 255

    img_B_intensive_G_uc = np.zeros([h, w])
    for node in B_intensive_G_uc:
        img_B_intensive_G_uc[node[0], node[1]] = 1
    G_uc_copy = G_uc.copy()

    G_uc1_copy = G_uc1.copy()
    for node in B_intensive_G_uc:
        nrb_node = []
        for edge in list(G_uc_copy.edges(nbunch=node)):
            nrb_node.append(edge[1])
        n = 0
        for node1 in nrb_node:
            if node1 not in B_intensive_G_uc:
                n += 1
        if n == len(nrb_node):
            continue
        (i, j) = node
        count, B_intensive_G_uc_sub = find_16x16(i, j, img_B_intensive_G_uc)
        if count > 10:
            count1, B_intensive_G_uc1_sub = find_16x16(i, j, image_zeros_copy / 255)
            for node_b in B_intensive_G_uc1_sub:
                if img_sk205[node_b[0], node_b[1]] != 255:
                    if node_b in G_uc1:
                        if G_uc1.degree[node_b] == 2:
                            nrb_node_b = []
                            for edge in list(G_uc1.edges(nbunch=node_b)):
                                nrb_node_b.append(edge[1])
                            G_uc1_copy.add_edge(nrb_node_b[0], nrb_node_b[-1])
                        G_uc1_copy.remove_node(node_b)

                        G_uc1.remove_node(node_b)


    image_zeros_new = np.zeros([h, w])
    for node in G_uc1.nodes:
        image_zeros_new[node[0], node[1]] = 255

    image_zeros_new_copy = np.zeros([h, w])
    for node in G_uc1_copy.nodes:
        image_zeros_new_copy[node[0], node[1]] = 255

    node_remove = []
    for node in G_uc_copy:
        if node not in G_uc1:
            node_remove.append(node)


    G_uc1 = optimize_2th(G_uc1=G_uc1, image_zeros_new=image_zeros_new, D_t=D_t, h=h)

    G_uc1_copy = optimize_2th(G_uc1=G_uc1_copy, image_zeros_new=image_zeros_new_copy, D_t=D_t, h=h)

    i1 = 0
    number = components.number_connected_components(G_uc1)
    for c in sorted(nx.connected_components(G_uc1), key=len, reverse=True):
        c = G_uc1.subgraph(c).copy()
        i1 += 1
        if i1 > (1 / ((int(number / 5)))) * number:
            score_veseel = 0
            for node0 in c.nodes:
                if (image_label[node0[0], node0[1], :] != [0, 0, 0]).any():
                    score_veseel += 1
            if score_veseel / len(c.nodes) < 0.3:
                for node in c.nodes:
                    G_uc1.remove_node(node)
                    if node in G_uc1_copy.nodes:
                        G_uc1_copy.remove_node(node)

    return G_uc1, D_t, image_zeros_205, node_remove, B_intensive_G_uc, G_uc_copy, G_uc1_copy