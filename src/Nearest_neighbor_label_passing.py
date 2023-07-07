import networkx as nx
import math

def cal_angle(point_a, point_b, point_c):
    a_x, b_x, c_x = point_a[0], point_b[0], point_c[0]
    a_y, b_y, c_y = point_a[1], point_b[1], point_c[1]
    if len(point_a) == len(point_b) == len(point_c) == 3:
        a_z, b_z, c_z = point_a[2], point_b[2], point_c[2]
    else:
        a_z, b_z, c_z = 0, 0, 0
    x1, y1, z1 = (a_x - b_x), (a_y - b_y), (a_z - b_z)
    x2, y2, z2 = (c_x - b_x), (c_y - b_y), (c_z - b_z)
    cos_b = (x1 * x2 + y1 * y2 + z1 * z2) / (
                math.sqrt(x1 ** 2 + y1 ** 2 + z1 ** 2) * (math.sqrt(x2 ** 2 + y2 ** 2 + z2 ** 2)) + 1e-8)
    B = math.degrees(math.acos(cos_b))
    return B

def segs_angle_othertwosege(segA, segB):
    angle1 = cal_angle(segA[-1], segA[0], segB[-1])
    angle2 = cal_angle(segA[-1], segA[1], segB[-1])
    angle3 = cal_angle(segA[-1], segB[1], segB[-1])
    angle = 0.8 * angle1 + 0.1 * angle2 + 0.1 * angle3
    return [angle, segA[1]]

def Find_Recent(x, y, x_, y_, k):
    list_stack_temp = []
    for i in range(len(x)):
        list_temp = []
        if x[i] != x_:
            length = math.sqrt(pow(x[i] - x_, 2) + pow(y[i] - y_, 2))
            if len(list_stack_temp) < k:
                list_stack_temp.append([(x[i], y[i]), length])
            else:
                for m in list_stack_temp:
                    list_temp.append(m[1])
                list_temp.append(length)
                list_temp.sort()
                if length != list_temp[-1]:
                    last_ = list_temp[-1]
                    for n in list_stack_temp:
                        if n[1] == last_:
                            list_stack_temp.remove(n)
                        else:
                            continue
                    list_stack_temp.append([(x[i], y[i]), length])
                else:
                    continue
        else:
            continue
    return list_stack_temp

def neighborhood(i, j, image, h, g, k1=0, k2=3):
    count = 0
    if i == 0:
        k1, k2 = 1, 3
    if i == h - 1:
        k1, k2 = 0, 2
    for hang in range(k1, k2):
        for lie in range(0, 3):
            if hang == lie == 1:
                continue
            if (image[i - 1 + hang, j - 1 + lie] != [0, 0, 0]).any():
                count += 1
                if i - 1 + hang < i:
                    continue
                elif i - 1 + hang == i and j - 1 + lie < j:
                    continue
                else:
                    g.add_edge((i, j), (i - 1 + hang, j - 1 + lie))
    return count, g


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
            if (image[i - 1 + hang, j - 1 + lie] != [0, 0, 0]).any():
                count += 1
                neighborhood_8.append((i - 1 + hang, j - 1 + lie))
    return count, neighborhood_8

def find_seg_passlabel(branchnode, nbr, graph):
    seg = []
    seg.append(branchnode)
    seg.append(nbr)
    if graph.degree[seg[-1]] != 2:
        return seg
    else:
        edge = [0, 0]
        while 1:
            b = seg[-1]
            for edge in list(graph.edges(nbunch=b)):
                if edge[1] not in seg:
                    seg.append(edge[1])
            if graph.degree[seg[-1]] != 2 or seg[-1] == b:
                break
        return seg

def nearest_neighbor_label_passing(h, w, image, G_direct_twoclasses, G_direct_v, G_direct_a, Disk_node, g4):
    g_vessel = nx.Graph()
    for i in range(0, h):
        for j in range(0, w):
            if (image[i, j, :] != [0, 0, 0]).any():
                g_vessel.add_node((i, j))
    for node in g_vessel.nodes:
        (i, j) = node
        k1, k2 = 0, 3
        count, g_vessel = neighborhood(i=i, j=j, image=image, h=h, g=g_vessel, k1=k1, k2=k2)

    g = nx.Graph()
    for node in G_direct_twoclasses:
        if G_direct_twoclasses.degree[node] > 2:
            in_degree = len(list(G_direct_twoclasses.in_edges(nbunch=node)))
            out_degree = len(list(G_direct_twoclasses.out_edges(nbunch=node)))
            if in_degree == 2 and out_degree == 1:
                in_nrb = []
                for edge in list(G_direct_twoclasses.in_edges(nbunch=node)):
                    in_nrb.append(edge[0])
                nrb_out = list(G_direct_twoclasses.out_edges(nbunch=node))[0][1]
                seg_out = find_seg_passlabel(branchnode=node, nbr=nrb_out, graph=g4)
                BN_next = seg_out[-1]
                in_degree_next = len(list(G_direct_twoclasses.in_edges(nbunch=BN_next)))
                out_degree_next = len(list(G_direct_twoclasses.out_edges(nbunch=BN_next)))
                nrb_out_next = []
                for edge in list(G_direct_twoclasses.out_edges(nbunch=BN_next)):
                    nrb_out_next.append(edge[1])
                j1 = 0
                for node_in in in_nrb:
                    seg_in = find_seg_passlabel(branchnode=node, nbr=node_in, graph=g4)
                    for i in range(2, len(seg_in)):
                        node_seg = seg_in[i]
                        if node_in in G_direct_a:
                            if (image[node_seg[0], node_seg[1], :] == [255, 0, 0]).all():
                                if i > 3 and i < 6:
                                    j1 += 1
                                    break
                        if node_in in G_direct_v:
                            if (image[node_seg[0], node_seg[1], :] == [0, 0, 255]).all():
                                if i > 3 and i < 6:
                                    j1 += 1
                                    break
                j2 = 0
                for node_out in nrb_out_next:
                    seg_out_next = find_seg_passlabel(branchnode=BN_next, nbr=node_out, graph=g4)
                    for i in range(2, len(seg_out_next)):
                        node_seg = seg_out_next[i]
                        if node_out in G_direct_a:
                            if (image[node_seg[0], node_seg[1], :] == [255, 0, 0]).all():
                                if i > 3 and i < 6:
                                    j2 += 1
                                    break
                        if node_out in G_direct_v:
                            if (image[node_seg[0], node_seg[1], :] == [0, 0, 255]).all():
                                if i > 3 and i < 6:
                                    j2 += 1
                                    break
                if in_degree_next == 1 and out_degree_next == 2:
                    i1 = 0
                    i2 = 0
                    for node_nrb in nrb_out_next:
                        if node_nrb in G_direct_a:
                            i1 += 1
                        if node_nrb in G_direct_v:
                            i2 += 1
                    if i1 == 1 and i2 == 1 and len(seg_out) < 12 * 2.5:
                        seg_out_mid = (
                        int((seg_out[0][0] + seg_out[-1][0]) / 2), int((seg_out[0][1] + seg_out[-1][1]) / 2))
                        l = math.sqrt((seg_out_mid[0] - seg_out[0][0]) ** 2 + (seg_out_mid[1] - seg_out[0][1]) ** 2)
                        if j1 == 0 and j2 == 0:
                            dist = math.sqrt(l ** 2 + 8 ** 2)
                            threshold = int(l) + 8
                        else:
                            dist = math.sqrt(l ** 2 + 4 ** 2)
                            threshold = int(l)
                        for i_ in range(-threshold, threshold):
                            for j_ in range(-threshold, threshold):
                                node_ = (seg_out_mid[0] + i_, seg_out_mid[1] + j_)
                                dist1 = math.sqrt((seg_out_mid[0] - node_[0]) ** 2 + (seg_out_mid[1] - node_[1]) ** 2)
                                if dist1 < dist:
                                    if (image[node_[0], node_[1], :] != [0, 0, 0]).any() or node_ == seg_out[
                                        0] or node_ == seg_out[-1]:
                                        g.add_node((node_[0], node_[1]))

            if in_degree + out_degree > 2 and node in G_direct_v and node in G_direct_a:
                in_nrb = []
                for edge in list(G_direct_twoclasses.in_edges(nbunch=node)):
                    in_nrb.append(edge[0])
                out_nrb = []
                for edge in list(G_direct_twoclasses.out_edges(nbunch=node)):
                    out_nrb.append(edge[1])
                j1 = 0
                for node_in in in_nrb:
                    seg_in = find_seg_passlabel(branchnode=node, nbr=node_in, graph=g4)
                    for i in range(2, len(seg_in)):
                        node_seg = seg_in[i]
                        if node_in in G_direct_a:
                            if (image[node_seg[0], node_seg[1], :] == [255, 0, 0]).all():
                                if i > 3 and i < 6:
                                    j1 += 1
                                    break
                        if node_in in G_direct_v:
                            if (image[node_seg[0], node_seg[1], :] == [0, 0, 255]).all():
                                if i > 3 and i < 6:
                                    j1 += 1
                                    break
                j2 = 0
                for node_out in out_nrb:
                    seg_out = find_seg_passlabel(branchnode=node, nbr=node_out, graph=g4)
                    for i in range(2, len(seg_out)):
                        node_seg = seg_out[i]
                        if node_out in G_direct_a:
                            if (image[node_seg[0], node_seg[1], :] == [255, 0, 0]).all():
                                if i > 3 and i < 6:
                                    j2 += 1
                                    break
                        if node_out in G_direct_v:
                            if (image[node_seg[0], node_seg[1], :] == [0, 0, 255]).all():
                                if i > 3 and i < 6:
                                    j2 += 1
                                    break
                if j1 == 0 and j2 == 0:
                    for i_1 in range(-8, 8):
                        for j_1 in range(-8, 8):
                            if (image[node[0] + i_1, node[1] + j_1, :] != [0, 0, 0]).any() and (node[0] + i_1, node[1] + j_1) not in g:
                                g.add_node((node[0] + i_1, node[1] + j_1))

            if in_degree == 2 and out_degree == 1 and node not in g:
                out_nrb = list(G_direct_twoclasses.out_edges(nbunch=node))[0][1]
                in_nrb = []
                for edge in list(G_direct_twoclasses.in_edges(nbunch=node)):
                    in_nrb.append(edge[0])
                in_seg1 = find_seg_passlabel(branchnode=node, nbr=in_nrb[0], graph=g4)
                in_seg2 = find_seg_passlabel(branchnode=node, nbr=in_nrb[1], graph=g4)
                (angle_in, nrb_) = segs_angle_othertwosege(in_seg1, in_seg2)
                if out_nrb in G_direct_a:
                    for node0 in in_nrb:
                        if node0 in G_direct_v:
                            seg = find_seg_passlabel(branchnode=node, nbr=node0, graph=g4)
                            remove_node = []
                            for node1 in seg:
                                if (image[node1[0], node1[1], :] == [0, 0, 255]).all():
                                    remove_node.append(node1)
                                if (image[node1[0], node1[1], :] != [0, 0, 255]).any():
                                    (i, j) = node1
                                    k1, k2 = 0, 3
                                    count, neighborhood_8 = Neighborhood_8(i=i, j=j, image=image, h=h, k1=k1, k2=k2)
                                    k0 = 0
                                    for node2 in neighborhood_8:
                                        if (image[node2[0], node2[1], :] == [0, 0, 255]).all():
                                            k0 += 1
                                    if k0 == 0:
                                        break
                                    else:
                                        remove_node.append(node1)
                            if len(remove_node) == 0:
                                continue
                            if remove_node[0] in G_direct_v:
                                G_direct_v.remove_node(remove_node[0])
                            if len(remove_node) > 3 and angle_in > 30:
                                continue
                            for i in range(1, len(remove_node)):
                                if remove_node[i] in G_direct_v:
                                    G_direct_v.remove_node(remove_node[i])
                                if i == 1:
                                    G_direct_a.add_edge(remove_node[i], remove_node[0])
                                else:
                                    G_direct_a.add_edge(remove_node[i], remove_node[i - 1])
                if out_nrb in G_direct_v:
                    for node0 in in_nrb:
                        if node0 in G_direct_a:
                            seg = find_seg_passlabel(branchnode=node, nbr=node0, graph=g4)
                            remove_node = []
                            for node1 in seg:
                                if (image[node1[0], node1[1], :] == [255, 0, 0]).all():
                                    remove_node.append(node1)
                                if (image[node1[0], node1[1], :] != [255, 0, 0]).any():
                                    (i, j) = node1
                                    k1, k2 = 0, 3
                                    count, neighborhood_8 = Neighborhood_8(i=i, j=j, image=image, h=h, k1=k1, k2=k2)
                                    k0 = 0
                                    for node2 in neighborhood_8:
                                        if (image[node2[0], node2[1], :] == [255, 0, 0]).all():
                                            k0 += 1
                                    if k0 == 0:
                                        break
                                    else:
                                        remove_node.append(node1)
                            if len(remove_node) == 0:
                                continue
                            if remove_node[0] in G_direct_a:
                                G_direct_a.remove_node(remove_node[0])
                            if len(remove_node) > 3 and angle_in > 30:
                                continue
                            for i in range(1, len(remove_node)):
                                G_direct_a.remove_node(remove_node[i])
                                if i == 1:
                                    G_direct_v.add_edge(remove_node[i], remove_node[0])
                                else:
                                    G_direct_v.add_edge(remove_node[i], remove_node[i - 1])
            if in_degree == 2 and out_degree == 2 and node not in g:
                in_nrb = []
                for edge in list(G_direct_twoclasses.in_edges(nbunch=node)):
                    in_nrb.append(edge[0])
                out_nrb = []
                for edge in list(G_direct_twoclasses.out_edges(nbunch=node)):
                    out_nrb.append(edge[1])
                out_seg1 = find_seg_passlabel(branchnode=node, nbr=out_nrb[0], graph=g4)
                out_seg2 = find_seg_passlabel(branchnode=node, nbr=out_nrb[1], graph=g4)
                (angle_out, nrb_) = segs_angle_othertwosege(out_seg1, out_seg2)
                in_seg1 = find_seg_passlabel(branchnode=node, nbr=in_nrb[0], graph=g4)
                in_seg2 = find_seg_passlabel(branchnode=node, nbr=in_nrb[1], graph=g4)
                (angle_in, nrb_) = segs_angle_othertwosege(in_seg1, in_seg2)
                if (image[node[0], node[1], :] == [255, 0, 0]).all():
                    for node0 in in_nrb:
                        if node0 in G_direct_a:
                            seg = find_seg_passlabel(branchnode=node, nbr=node0, graph=g4)
                            remove_node = []
                            for node1 in seg:
                                if (image[node1[0], node1[1], :] == [255, 0, 0]).all():
                                    remove_node.append(node1)
                                if (image[node1[0], node1[1], :] != [255, 0, 0]).any():
                                    (i, j) = node1
                                    k1, k2 = 0, 3
                                    count, neighborhood_8 = Neighborhood_8(i=i, j=j, image=image, h=h, k1=k1, k2=k2)
                                    k0 = 0
                                    for node2 in neighborhood_8:
                                        if (image[node2[0], node2[1], :] == [255, 0, 0]).all():
                                            k0 += 1
                                    if k0 == 0:
                                        break
                                    else:
                                        remove_node.append(node1)
                            if remove_node[0] in G_direct_a:
                                G_direct_a.remove_node(remove_node[0])
                            if len(remove_node) > 3 and angle_in > 30:
                                continue
                            for i in range(1, len(remove_node)):
                                G_direct_a.remove_node(remove_node[i])
                                if i == 1:
                                    G_direct_v.add_edge(remove_node[i], remove_node[0])
                                else:
                                    G_direct_v.add_edge(remove_node[i], remove_node[i - 1])
                    for node0 in out_nrb:
                        if node0 in G_direct_a:
                            seg = find_seg_passlabel(branchnode=node, nbr=node0, graph=g4)
                            remove_node = []
                            for node1 in seg:
                                if (image[node1[0], node1[1], :] == [255, 0, 0]).all():
                                    remove_node.append(node1)
                                if (image[node1[0], node1[1], :] != [255, 0, 0]).any():
                                    (i, j) = node1
                                    k1, k2 = 0, 3
                                    count, neighborhood_8 = Neighborhood_8(i=i, j=j, image=image, h=h, k1=k1, k2=k2)
                                    k0 = 0
                                    for node2 in neighborhood_8:
                                        if (image[node2[0], node2[1], :] == [255, 0, 0]).all():
                                            k0 += 1
                                    if k0 == 0:
                                        break
                                    else:
                                        remove_node.append(node1)
                            if remove_node[0] in G_direct_a:
                                G_direct_a.remove_node(remove_node[0])
                            if len(remove_node) > 3 and angle_out > 30:
                                continue
                            for i in range(1, len(remove_node)):
                                G_direct_a.remove_node(remove_node[i])
                                if i == 1:
                                    G_direct_v.add_edge(remove_node[i], remove_node[0])
                                else:
                                    G_direct_v.add_edge(remove_node[i], remove_node[i - 1])
                if (image[node[0], node[1], :] == [0, 0, 255]).all():
                    for node0 in in_nrb:
                        if node0 in G_direct_v:
                            seg = find_seg_passlabel(branchnode=node, nbr=node0, graph=g4)
                            remove_node = []
                            for node1 in seg:
                                if (image[node1[0], node1[1], :] == [0, 0, 255]).all():
                                    remove_node.append(node1)
                                if (image[node1[0], node1[1], :] != [0, 0, 255]).any():
                                    (i, j) = node1
                                    k1, k2 = 0, 3
                                    count, neighborhood_8 = Neighborhood_8(i=i, j=j, image=image, h=h, k1=k1, k2=k2)
                                    k0 = 0
                                    for node2 in neighborhood_8:
                                        if (image[node2[0], node2[1], :] == [0, 0, 255]).all():
                                            k0 += 1
                                    if k0 == 0:
                                        break
                                    else:
                                        remove_node.append(node1)
                            if remove_node[0] in G_direct_v:
                                G_direct_v.remove_node(remove_node[0])
                            if len(remove_node) > 3 and angle_in > 30:
                                continue
                            for i in range(1, len(remove_node)):
                                if remove_node[i] not in G_direct_v.nodes:
                                    continue
                                G_direct_v.remove_node(remove_node[i])
                                if i == 1:
                                    G_direct_a.add_edge(remove_node[i], remove_node[0])
                                else:
                                    G_direct_a.add_edge(remove_node[i], remove_node[i - 1])
                    for node0 in out_nrb:
                        if node0 in G_direct_v:
                            seg = find_seg_passlabel(branchnode=node, nbr=node0, graph=g4)
                            remove_node = []
                            for node1 in seg:
                                if (image[node1[0], node1[1], :] == [0, 0, 255]).all():
                                    remove_node.append(node1)
                                if (image[node1[0], node1[1], :] != [0, 0, 255]).any():
                                    (i, j) = node1
                                    k1, k2 = 0, 3
                                    count, neighborhood_8 = Neighborhood_8(i=i, j=j, image=image, h=h, k1=k1, k2=k2)
                                    k0 = 0
                                    for node2 in neighborhood_8:
                                        if (image[node2[0], node2[1], :] == [0, 0, 255]).all():
                                            k0 += 1
                                    if k0 == 0:
                                        break
                                    else:
                                        remove_node.append(node1)
                            if remove_node[0] in G_direct_a:
                                G_direct_a.remove_node(remove_node[0])
                            if len(
                                    remove_node) > 3 and angle_out > 30:
                                continue
                            for i in range(1, len(remove_node)):
                                if remove_node[i] not in G_direct_a.nodes:
                                    continue
                                G_direct_a.remove_node(remove_node[i])
                                if i == 1:
                                    G_direct_v.add_edge(remove_node[i], remove_node[0])
                                else:
                                    G_direct_v.add_edge(remove_node[i], remove_node[i - 1])

    B = G_direct_twoclasses
    while 1:
        label_set = B
        if len(B) == 0:
            break
        set_next = []
        for node in label_set:
            (i, j) = node
            k1, k2 = 0, 3
            count, neighborhood_8 = Neighborhood_8(i=i, j=j, image=image, h=h, k1=k1, k2=k2)
            for node1 in neighborhood_8:
                if node not in g or node1 not in g:
                    if node in G_direct_a and node1 not in G_direct_v:
                        if (image[node[0], node[1], :] != [0, 0,
                                                           255]).any() and node in g_vessel and node not in g:
                            image[node[0], node[1], :] = [0, 0, 255]
                        if node1 in G_direct_a or node1 in Disk_node or node1 in g:
                            continue
                        image[node1[0], node1[1], :] = [0, 0, 255]
                        set_next.append(node1)
                        G_direct_a.add_node(node1)
                    if node in G_direct_v and node1 not in G_direct_a:
                        if (image[node[0], node[1], :] != [255, 0,
                                                           0]).any() and node in g_vessel and node not in g:
                            image[node[0], node[1], :] = [255, 0, 0]
                        if node1 in G_direct_v or node1 in Disk_node or node1 in g:
                            continue
                        image[node1[0], node1[1], :] = [255, 0, 0]
                        G_direct_v.add_node(node1)
                        set_next.append(node1)
        B = set_next
    return image