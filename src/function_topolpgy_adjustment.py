import math

def find_nbrBN(graph, node, next_nbr_bn, next_nbr_norm, Bn_g3):
    for edge in list(graph.edges(nbunch=node)):
        if edge[1] in Bn_g3 and edge[1] not in next_nbr_bn:
            next_nbr_bn.append(edge[1])
        if graph.degree[edge[1]] < 3:
            next_nbr_norm.append(edge[1])
    return next_nbr_bn, next_nbr_norm

def find_nbrBN_new(graph, node, next_nbr_bn, next_nbr_norm, Bn_g3):
    for edge in list(graph.edges(nbunch=node)):
        if edge[1] in Bn_g3 and edge[1] not in next_nbr_bn:
            next_nbr_bn.append(edge[1])
        if graph.degree[edge[1]] < 3 and edge[1] not in next_nbr_norm:
            next_nbr_norm.append(edge[1])
    return next_nbr_bn, next_nbr_norm

def find_BN_NN_BN(graph, node, bn_nn_bn, addedge_node):
    bn_new = []
    for edge in list(graph.edges(nbunch=node)):
        if graph.degree[edge[1]] == 2:
            for edge1 in list(graph.edges(nbunch=edge[1])):
                if edge1[1] != edge[0]:
                    if graph.degree[edge1[1]] > 2:
                        bn_nn_bn.append(edge1[1])
                        bn_new.append(edge[1])
                        for edge2 in list(graph.edges(nbunch=edge1[1])):
                            if edge2[1] != bn_new[-1]:
                                addedge_node.append(edge2[1])
        addedge_node.append(edge[1])
    addedge_node_new = []
    if len(bn_new) != 0:
        for node in addedge_node:
            if node != bn_new[-1]:
                addedge_node_new.append(node)
    return bn_nn_bn, addedge_node_new, bn_new

def find_BN_NN_BN_new(graph, node, bn_nn_bn, addedge_node):
    bn_new = []
    for edge in list(graph.edges(nbunch=node)):
        if graph.degree[edge[1]] == 2:
            for edge1 in list(graph.edges(nbunch=edge[1])):
                if edge1[1] != edge[0]:
                    if graph.degree[edge1[1]] == 3:
                        bn_nn_bn.append(edge1[1])
                        bn_new.append(edge[1])
                        for edge2 in list(graph.edges(nbunch=edge1[1])):
                            if edge2[1] != bn_new[-1] and edge2[1] not in addedge_node:
                                addedge_node.append(edge2[1])
        addedge_node.append(edge[1])
    addedge_node_new = []
    if len(bn_new) != 0:
        for node in addedge_node:
            if node != bn_new[-1]:
                addedge_node_new.append(node)
    return bn_nn_bn, addedge_node_new, bn_new

def find_seg(branchnode, nbr, graph):
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
            if graph.degree[seg[-1]] != 2:
                break
        return seg

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

def topology_adjustment1(g3):
    Bn_g3 = []
    g4 = g3.copy()
    for node in g3.nodes:
        if g3.degree[node] > 2:
            Bn_g3.append(node)
    while 1:
        next_nbr_bn = []
        next_nbr_norm = []
        if len(Bn_g3) == 0:
            break
        node = Bn_g3.pop(-1)
        if node not in g4:
            continue
        next_nbr_bn.append(node)
        while 1:
            next_nbr_bn, next_nbr_norm = find_nbrBN(g4, node, next_nbr_bn, next_nbr_norm, Bn_g3)
            n1 = len(next_nbr_bn)
            for node in next_nbr_bn:
                if node == next_nbr_bn[0]:
                    continue
                next_nbr_bn, next_nbr_norm = find_nbrBN(g4, node, next_nbr_bn, next_nbr_norm, Bn_g3)
            n2 = len(next_nbr_bn)
            if n1 == n2:
                break
            for i in range(n1, len(next_nbr_bn)):
                node = next_nbr_bn[i]
                next_nbr_bn, next_nbr_norm = find_nbrBN(g4, node, next_nbr_bn, next_nbr_norm, Bn_g3)
            if n2 == len(next_nbr_bn):
                break
        if len(next_nbr_bn) != 1:
            node_1 = 0
            node_2 = 0
            for i in range(0, len(next_nbr_bn)):
                node_1 += next_nbr_bn[i][0]
                node_2 += next_nbr_bn[i][1]
                g4.remove_node(next_nbr_bn[i])
            node_new = (int(node_1 / len(next_nbr_bn)), int(node_2 / len(next_nbr_bn)))
            g4.add_node(node_new)
            for node in next_nbr_norm:
                g4.add_edge(node, node_new)
        if len(Bn_g3) == 0:
            break
    return g4

def topology_adjustment2(g4):
    Bn_g4 = []
    for node in g4.nodes:
        if g4.degree[node] == 3:
            Bn_g4.append(node)
    while 1:
        bn_nn_bn = []
        addedge_node = []
        if len(Bn_g4) == 0:
            break
        node = Bn_g4.pop(-1)
        if node not in g4:
            continue
        bn_nn_bn.append(node)
        bn_nn_bn, addedge_node_1, bn_new_1 = find_BN_NN_BN_new(g4, node, bn_nn_bn, addedge_node)
        if len(bn_nn_bn) == 2:
            bn_nn_bn_next = [bn_nn_bn[-1]]
            addedge_node = []
            bn_nn_bn_next, addedge_node_2, bn_new_2 = find_BN_NN_BN_new(g4, bn_nn_bn[-1], bn_nn_bn_next, addedge_node)
            bn_nn_bn_next_copy = bn_nn_bn_next.copy()
            if len(bn_nn_bn_next) == 3:
                bn_h = [bn_nn_bn_next[0][0], bn_nn_bn_next[1][0], bn_nn_bn_next[2][0]]
                if min(bn_h) == bn_nn_bn_next[1][0]:
                    bn_nn_bn_new = [bn_nn_bn_next[0], bn_nn_bn_next[1]]
                    bn_new = [bn_new_2[0]]
                    addedge_node = []
                    bn_nn_bn_next = [bn_nn_bn_new[-1]]
                    bn_nn_bn_next_, addedge_node_new, bn_new_ = find_BN_NN_BN_new(g4, bn_nn_bn_new[-1], bn_nn_bn_next,
                                                                                addedge_node)
                elif min(bn_h) == bn_nn_bn_next_copy[2][0]:
                    bn_nn_bn_new = [bn_nn_bn_next[0], bn_nn_bn_next[2]]
                    bn_new = [bn_new_2[1]]
                    addedge_node = []
                    bn_nn_bn_next = [bn_nn_bn_new[-1]]
                    bn_nn_bn_next_, addedge_node_new, bn_new_ = find_BN_NN_BN_new(g4, bn_nn_bn_new[-1], bn_nn_bn_next,
                                                                                  addedge_node)
                else:
                    bn_w = [bn_nn_bn_next[0][1], bn_nn_bn_next[1][1], bn_nn_bn_next[2][1]]
                    if min(bn_w) == bn_nn_bn_next[1][1]:
                        bn_nn_bn_new = [bn_nn_bn_next[0], bn_nn_bn_next[1]]
                        bn_new = [bn_new_2[0]]
                        addedge_node = []
                        bn_nn_bn_next = [bn_nn_bn_new[-1]]
                        bn_nn_bn_next_, addedge_node_new, bn_new_ = find_BN_NN_BN_new(g4, bn_nn_bn_new[-1], bn_nn_bn_next,
                                                                                      addedge_node)
                    if min(bn_w) == bn_nn_bn_next_copy[2][1]:
                        bn_nn_bn_new = [bn_nn_bn_next[0], bn_nn_bn_next[2]]
                        bn_new = [bn_new_2[1]]
                        addedge_node = []
                        bn_nn_bn_next = [bn_nn_bn_new[-1]]
                        bn_nn_bn_next_, addedge_node_new, bn_new_ = find_BN_NN_BN_new(g4, bn_nn_bn_new[-1], bn_nn_bn_next,
                                                                                      addedge_node)
            else:
                bn_nn_bn_new = bn_nn_bn
                bn_new = bn_new_1
                addedge_node_new = addedge_node_1

            g4.remove_node(bn_nn_bn_new[0])
            g4.remove_node(bn_nn_bn_new[1])
            for node in addedge_node_new:
                g4.add_edge(bn_new[0], node)
        if len(Bn_g4) == 0:
            break
    return g4


def topology_adjustment3(g3):
    Bn_g3 = []
    g4 = g3.copy()
    for node in g3.nodes:
        if g3.degree[node] > 2:
            Bn_g3.append(node)
    last_nbr_bn = []
    while 1:
        next_nbr_bn = []
        next_nbr_norm = []
        if len(Bn_g3) == 0:
            break
        node = Bn_g3.pop(-1)
        if node in last_nbr_bn:
            continue
        if node not in g4:
            continue
        next_nbr_bn.append(node)
        while 1:
            next_nbr_bn, next_nbr_norm = find_nbrBN_new(g4, node, next_nbr_bn, next_nbr_norm, Bn_g3)
            n1 = len(next_nbr_bn)
            for node in next_nbr_bn:
                if node == next_nbr_bn[0]:
                    continue
                next_nbr_bn, next_nbr_norm = find_nbrBN_new(g4, node, next_nbr_bn, next_nbr_norm, Bn_g3)
            n2 = len(next_nbr_bn)
            if n1 == n2:
                break
            for i in range(n1, len(next_nbr_bn)):
                node = next_nbr_bn[i]
                next_nbr_bn, next_nbr_norm = find_nbrBN_new(g4, node, next_nbr_bn, next_nbr_norm, Bn_g3)
            if n2 == len(next_nbr_bn):
                break
        for node1 in next_nbr_bn:
            last_nbr_bn.append(node1)
        if len(next_nbr_bn) == 2:
            if len(next_nbr_norm) == 3:
                for node1 in next_nbr_norm:
                    nrb = []
                    for edge in list(g4.edges(nbunch=node1)):
                        nrb.append(edge[1])
                    k = 0
                    for node0 in nrb:
                        if node0 in next_nbr_bn:
                            k += 1
                    if k == len(nrb):
                        next_nrb_two = []
                        for node in next_nbr_norm:
                            if node != node1:
                                next_nrb_two.append(node)
                        break
                two_seg = []
                for nrb0 in next_nrb_two:
                    for edge in list(g4.edges(nbunch=nrb0)):
                        if edge[1] not in next_nbr_bn:
                            continue
                        else:
                            bn = edge[1]
                    nrb0_seg = find_seg(branchnode=bn, nbr=nrb0, graph=g4)
                    two_seg.append(nrb0_seg)
                [angle, nrb_] = segs_angle_othertwosege(two_seg[0], two_seg[-1])
                if abs(angle-180) < 60:
                    continue
                node_1 = 0
                node_2 = 0
                for i in range(0, len(next_nbr_bn)):
                    node_1 += next_nbr_bn[i][0]
                    node_2 += next_nbr_bn[i][1]
                    g4.remove_node(next_nbr_bn[i])
                node_new = (int(node_1 / len(next_nbr_bn)), int(node_2 / len(next_nbr_bn)))
                g4.add_node(node_new)
                for node in next_nbr_norm:
                    g4.add_edge(node, node_new)

        if len(Bn_g3) == 0:
            break
    return g4