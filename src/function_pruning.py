import networkx as nx
import math

def find_seg_new(branchnode, nbr, close_disc, graph, Disk_node, B_g4):
    seg = []
    seg.append(branchnode)
    seg.append(nbr)
    if graph.degree[seg[-1]] != 2:
        return seg, close_disc
    else:
        edge = [0, 0]
        while 1:
            b = seg[-1]
            m = 0
            n = 0
            for edge in list(graph.edges(nbunch=b)):
                if edge[1] not in seg:
                    n += 1
                    if edge[1] not in Disk_node:
                        seg.append(edge[1])
                    if edge[1] in Disk_node:
                        m = 1
                        if edge[0] not in close_disc:
                            close_disc.append(edge[0])
            if m != 0:
                break
            if n == 0:
                break

            for edge in list(graph.edges(nbunch=seg[-1])):
                if edge[1] in Disk_node:
                    if edge[0] not in close_disc:
                        close_disc.append(edge[0])
                    break
                break
            if seg[-1] in B_g4 or edge[1] in Disk_node:
                break
        return seg, close_disc

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
        if x[i] != x_ or y[i] != y_:
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

def Pruning(g4, Disk_node, image_label, ratio = 1):
    D_t = []
    close_disc = []
    nrb_1_g4 = []
    nrb_bn_g4 = []
    for node in g4:
        if g4.degree[node] == 1:
            nrb_1_g4.append(node)
            D_t.append(node)
        if g4.degree[node] > 2:
            nrb_bn_g4.append(node)
    B_g4_copy = nrb_1_g4 + nrb_bn_g4
    dt_remove = []
    i = 0
    for node in D_t:
        i += 1
        if node in dt_remove:
            continue

        nrb = list(g4.edges(nbunch=node))[0][1]
        seg, close_disc = find_seg_new(node, nrb, close_disc, g4, Disk_node, B_g4_copy)

        degree = g4.degree[seg[-1]]
        if degree > 3 and len(seg) < 6 * ratio:
            for i in range(0, len(seg) - 1):
                g4.remove_node(seg[i])
        if degree == 3:
            nrb_seglast = []
            for edge in list(g4.edges(nbunch=seg[-1])):
                nrb_seglast.append(edge[1])
            seg_next = []
            for node1 in nrb_seglast:
                seg1, close_disc = find_seg_new(seg[-1], node1, close_disc, g4, Disk_node, B_g4_copy)
                if seg1[-1] != seg[0]:
                    seg_next.append(seg1)
            seg_1 = [seg_next[0][0], seg_next[0][1], seg_next[0][int(len(seg_next[0]) / 2)], seg_next[0][-2],
                     seg_next[0][-1]]
            seg_2 = [seg_next[1][0], seg_next[1][1], seg_next[1][int(len(seg_next[1]) / 2)], seg_next[1][-2],
                     seg_next[1][-1]]
            [angle, nrb0] = segs_angle_othertwosege(seg_1, seg_2)

            if len(seg) == 2 and (g4.degree[seg_1[-1]] == 1 and len(seg_next[0]) == 2):
                for i in range(0, len(seg) - 1):
                    g4.remove_node(seg[i])
                D_t = [x for x in D_t if x != seg_1[-1]]
                continue
            if len(seg) == 2 and (g4.degree[seg_2[-1]] == 1 and len(seg_next[1]) == 2):
                for i in range(0, len(seg) - 1):
                    g4.remove_node(seg[i])
                D_t = [x for x in D_t if x != seg_2[-1]]
                continue

            if angle < 120:
                continue
            dist = math.sqrt((seg[0][0] - seg[-1][0]) ** 2 + (seg[0][1] - seg[-1][1]) ** 2)
            if len(seg) < 6 * ratio and dist < 10 * 2.5:
                for i in range(0, len(seg) - 1):
                    g4.remove_node(seg[i])
                continue
            score_veseel = 0
            for node0 in seg:
                if (image_label[node0[0], node0[1], :] != [0, 0, 0]).any():
                    score_veseel += 1
            if score_veseel / len(seg) < 0.3:
                for i in range(0, len(seg) - 1):
                    g4.remove_node(seg[i])
                dt_remove.append(seg[-1])
                continue
        if degree == 1:
            score_veseel = 0
            for node0 in seg:
                if (image_label[node0[0], node0[1], :] != [0, 0, 0]).any():
                    score_veseel += 1
            if g4.degree[seg[-1]] == 1 and ((len(seg) > 12 * ratio and score_veseel / len(seg) < 0.3) or (
                    len(seg) < 12 * ratio and score_veseel / len(seg) < 0.5)):
                for i in range(0, len(seg) - 1):
                    g4.remove_node(seg[i])
                dt_remove.append(seg[-1])
                continue
    return g4

def find_seg_passlabel_new(branchnode, nbr, graph, close_disc):
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
            if graph.degree[seg[-1]] != 2 or seg[-1] in close_disc:
                break
        return seg

def Connecting(g4, Disk_node, K = 1):
    close_disc = []
    nrb_1_g4 = []
    nrb_bn_g4 = []
    D_t = []
    for node in g4:
        if g4.degree[node] == 1:
            nrb_1_g4.append(node)
            D_t.append(node)
        if g4.degree[node] > 2:
            nrb_bn_g4.append(node)
    B_g4_copy = nrb_1_g4 + nrb_bn_g4
    k = 0
    while 1:
        k += 1
        j = 0
        i = 0
        for c in sorted(nx.connected_components(g4), key=len, reverse=True):
            c = g4.subgraph(c).copy()
            i += 1
            if i == 1:
                g4_max = c
                X = []
                Y = []
                for node0 in g4_max:
                    X.append(node0[1])
                    Y.append(node0[0])
            num_g4_max = len(g4_max.nodes)
            if i != 1:
                if len(c.nodes) < 15:
                    continue
                dist0 = 1000
                node_ = (0, 0)
                node_dt = []
                for node in c:
                    if c.degree[node] == 1:
                        node_dt.append(node)
                for node in node_dt:
                    [[(x, y), dist]] = Find_Recent(x=X, y=Y, x_=node[1], y_=node[0], k=1)
                    if dist < dist0:
                        dist0 = dist
                        [node_, node_new] = [(y, x), node]
                if dist0 < 13 * 2.5:
                    j += 1
                    if g4.degree[node_] == 1:
                        seg1, close_disc = find_seg_new(node, list(g4.edges(nbunch=node))[0][1], close_disc, g4, Disk_node, B_g4_copy)
                        seg2, close_disc = find_seg_new(node_, list(g4.edges(nbunch=node_))[0][1], close_disc, g4, Disk_node, B_g4_copy)
                        [angle, nrb0] = segs_angle_othertwosege(seg1, seg2)
                        if angle > 120:
                            g4.add_edge(node_, node_new)
                    else:
                        g4.add_edge(node_, node_new)

        if k == K:
            break
    return g4

def Connecting0(g4, Disk_node, Disc_size, H, K = 1, K_s=1):
    close_disc = []
    nrb_1_g4 = []
    nrb_bn_g4 = []
    D_t = []
    for node in g4:
        if g4.degree[node] == 1:
            nrb_1_g4.append(node)
            D_t.append(node)
        if g4.degree[node] > 2:
            nrb_bn_g4.append(node)
    B_g4_copy = nrb_1_g4 + nrb_bn_g4
    k = 0
    while 1:
        k += 1
        j = 0
        i = 0
        for c in sorted(nx.connected_components(g4), key=len, reverse=True):
            c = g4.subgraph(c).copy()
            i += 1
            if i == 1:
                g4_max = c
                X = []
                Y = []
                for node0 in g4_max:
                    X.append(node0[1])
                    Y.append(node0[0])
            num_g4_max = len(g4_max.nodes)
            if i != 1:
                if len(c.nodes) < 15:
                    continue
                dist0 = 1000
                node_ = (0, 0)
                node_dt = []
                for node in c:
                    if c.degree[node] == 1:
                        node_dt.append(node)
                for node in node_dt:
                    [[(x, y), dist]] = Find_Recent(x=X, y=Y, x_=node[1], y_=node[0], k=1)
                    if dist < dist0:
                        dist0 = dist
                        [node_, node_new] = [(y, x), node]
                dist1 = math.sqrt((node_new[0] - H[0]) ** 2 + (node_new[1] - H[1]) ** 2)
                if dist1 > Disc_size * 3:
                    if dist0 < 13 * K_s:
                        j += 1
                        if g4.degree[node_] == 1:
                            seg1, close_disc = find_seg_new(node_new, list(g4.edges(nbunch=node_new))[0][1], close_disc, g4,
                                                            Disk_node, B_g4_copy)
                            seg2, close_disc = find_seg_new(node_, list(g4.edges(nbunch=node_))[0][1], close_disc, g4, Disk_node, B_g4_copy)
                            [angle, nrb0] = segs_angle_othertwosege(seg1, seg2)
                            if angle > 120:
                                g4.add_edge(node_, node_new)
                        else:
                            g4.add_edge(node_, node_new)
        if k == K:
            break
    return g4

def Connecting1(G_direct_new, g4, g4_not_cut, starting_point, Disc_size, H, close_disc, K = 1):
    k = 0
    while 1:
        k += 1
        j = 0
        for node in starting_point:
            G_direct_new_copy = G_direct_new.copy()
            dist1 = math.sqrt((node[0] - H[0]) ** 2 + (node[1] - H[1]) ** 2)
            if dist1 > Disc_size * 3 and dist1 < Disc_size * 12:
                seg = find_seg_passlabel_new(branchnode=node, nbr=list(G_direct_new.out_edges(nbunch=node))[0][1],
                                             graph=g4_not_cut, close_disc=close_disc)
                for node0 in seg:
                    G_direct_new_copy = [x for x in G_direct_new_copy if x != node0]
                X = []
                Y = []
                for node1 in G_direct_new_copy:
                    X.append(node1[1])
                    Y.append(node1[0])
                [[(x, y), dist]] = Find_Recent(x=X, y=Y, x_=node[1], y_=node[0], k=1)
                if dist < 13 * 2.5:
                    j += 1
                    G_direct_new.add_edge((y, x), node)
                    g4.add_edge((y, x), node)
                    g4_not_cut.add_edge((y, x), node)
                    starting_point = [x for x in starting_point if x != node]
                    nrb_digraph_out = []
                    for edge in list(G_direct_new.out_edges(nbunch=(y, x))):
                        if edge[1] != node:
                            nrb_digraph_out.append(edge[1])
                    nrb_digraph_in = list(G_direct_new.in_edges(nbunch=(y, x)))[0][0]
                    if len(nrb_digraph_out) == 0:
                        continue
                    else:
                        for edge in list(g4_not_cut.edges(nbunch=(y, x))):
                            seg_in_all = find_seg_passlabel_new(branchnode=(y, x), nbr=edge[1],
                                                                graph=g4_not_cut, close_disc=close_disc)
                            if nrb_digraph_out[0] in seg_in_all and seg_in_all[1] != nrb_digraph_out:
                                G_direct_new.remove_edge((y, x), nrb_digraph_out[0])
                                G_direct_new.add_edge((y, x), seg_in_all[1])
                                G_direct_new.add_edge(seg_in_all[1], nrb_digraph_out[0])
                            if nrb_digraph_in in seg_in_all and seg_in_all[1] != nrb_digraph_in:
                                G_direct_new.remove_edge(nrb_digraph_in, (y, x))
                                G_direct_new.add_edge(nrb_digraph_in, seg_in_all[1])
                                G_direct_new.add_edge(seg_in_all[1], (y, x))
        if k == K:
            break
    return G_direct_new, g4, g4_not_cut, starting_point

def Connecting2(G_direct_new, g4, g4_not_cut, starting_point, Disc_size, H, close_disc, image_label, K = 1):
    k = 0
    while 1:
        k += 1
        j = 0
        for node in starting_point:
            G_direct_new_copy = G_direct_new.copy()
            dist1 = math.sqrt((node[0] - H[0]) ** 2 + (node[1] - H[1]) ** 2)
            if dist1 > Disc_size * 3 and dist1 < Disc_size * 15:
                seg = find_seg_passlabel_new(branchnode=node, nbr=list(G_direct_new.out_edges(nbunch=node))[0][1],
                                             graph=g4_not_cut, close_disc=close_disc)
                for node0 in seg:
                    G_direct_new_copy = [x for x in G_direct_new_copy if x != node0]
                X = []
                Y = []
                for node1 in G_direct_new_copy:
                    X.append(node1[1])
                    Y.append(node1[0])
                [[(x, y), dist]] = Find_Recent(x=X, y=Y, x_=node[1], y_=node[0], k=1)
                if dist < 13 * 2.5:
                    nrb = []
                    for edge in list(g4_not_cut.edges(nbunch=(y, x))):
                        nrb.append(edge[1])
                    seg1 = find_seg_passlabel_new(branchnode=(y, x), nbr=nrb[0],
                                                  graph=g4_not_cut, close_disc=close_disc)
                    seg2 = find_seg_passlabel_new(branchnode=(y, x), nbr=nrb[1],
                                                  graph=g4_not_cut, close_disc=close_disc)
                    seg = find_seg_passlabel_new(branchnode=node, nbr=list(g4_not_cut.edges(nbunch=node))[0][1],
                                                  graph=g4_not_cut, close_disc=close_disc)
                    score_v1 = 0
                    score_a1 = 0
                    score_v2 = 0
                    score_a2 = 0
                    score_v = 0
                    score_a = 0
                    for node_ in seg1:
                        if (image_label[node_[0], node_[1], :] == [255, 0, 0]).all():
                            score_v1 += 1
                        if (image_label[node_[0], node_[1], :] == [0, 0, 255]).all():
                            score_a1 += 1
                    for node_ in seg2:
                        if (image_label[node_[0], node_[1], :] == [255, 0, 0]).all():
                            score_v2 += 1
                        if (image_label[node_[0], node_[1], :] == [0, 0, 255]).all():
                            score_a2 += 1
                    for node_ in seg:
                        if (image_label[node_[0], node_[1], :] == [255, 0, 0]).all():
                            score_v += 1
                        if (image_label[node_[0], node_[1], :] == [0, 0, 255]).all():
                            score_a += 1
                    if (score_a1 / (score_a1 + score_v1 + 1e-6) > score_v1 / (score_a1 + score_v1 + 1e-6)) and (
                            score_a2 / (score_a2 + score_v2 + 1e-6) > score_v2 / (score_a1 + score_v2 + 1e-6)) and \
                            (score_a / (score_a + score_v + 1e-6) < score_v / (score_a + score_v + 1e-6)):
                        continue
                    if (score_a1 / (score_a1 + score_v1 + 1e-6) < score_v1 / (score_a1 + score_v1 + 1e-6)) and (
                            score_a2 / (score_a2 + score_v2 + 1e-6) < score_v2 / (score_a1 + score_v2 + 1e-6)) and \
                            (score_a / (score_a + score_v + 1e-6) > score_v / (score_a + score_v + 1e-6)):
                        continue

                    j += 1
                    G_direct_new.add_edge((y, x), node)
                    g4.add_edge((y, x), node)
                    g4_not_cut.add_edge((y, x), node)
                    starting_point = [x for x in starting_point if x != node]
                    nrb_digraph_out = []
                    for edge in list(G_direct_new.out_edges(nbunch=(y, x))):
                        if edge[1] != node:
                            nrb_digraph_out.append(edge[1])
                    nrb_digraph_in = list(G_direct_new.in_edges(nbunch=(y, x)))[0][0]

                    if len(nrb_digraph_out) == 0:
                        continue
                    else:
                        for edge in list(g4_not_cut.edges(nbunch=(y, x))):
                            seg_in_all = find_seg_passlabel_new(branchnode=(y, x), nbr=edge[1],
                                                                graph=g4_not_cut, close_disc=close_disc)
                            if nrb_digraph_out[0] in seg_in_all and seg_in_all[1] != nrb_digraph_out:
                                G_direct_new.remove_edge((y, x), nrb_digraph_out[0])
                                G_direct_new.add_edge((y, x), seg_in_all[1])
                                G_direct_new.add_edge(seg_in_all[1], nrb_digraph_out[0])
                            if nrb_digraph_in in seg_in_all and seg_in_all[1] != nrb_digraph_in:
                                G_direct_new.remove_edge(nrb_digraph_in, (y, x))
                                G_direct_new.add_edge(nrb_digraph_in, seg_in_all[1])
                                G_direct_new.add_edge(seg_in_all[1], (y, x))
        if k == K:
            break
    return G_direct_new, g4, g4_not_cut, starting_point