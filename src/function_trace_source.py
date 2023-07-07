import networkx as nx
import math

def find_seg(branchnode, nbr, graph):
    seg = []
    seg.append(branchnode)
    seg.append(nbr)
    if graph.degree[seg[-1]] != 2:
        return seg
    else:
        while 1:
            b = seg[-1]
            for edge in list(graph.edges(nbunch=b)):
                if edge[1] not in seg:
                    seg.append(edge[1])
            if graph.degree[seg[-1]] != 2:
                break
        return seg

def find_10x10(i, j, image):
    score_v = 0
    score_a = 0
    for m in range(-5, 5):
        for n in range(-5, 5):
            if (image[i + m, j + n, :] == [255, 0, 0]).all():
                score_v += 1
            if (image[i + m, j + n, :] == [0, 0, 255]).all():
                score_a += 1
    return score_v, score_a

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

def segs_angle(segA, segB):
    angle1 = cal_angle(segA[-1], segA[0], segB[0])
    angle2 = cal_angle(segA[-1], segA[1], segB[0])
    angle3 = cal_angle(segA[-1], segB[-2], segB[0])
    angle = 0.8 * angle1 + 0.1 * angle2 + 0.1 * angle3
    return [angle, segA[1]]

def segs_angle_new(segA, segB, graph):
    angle1 = cal_angle(segA[-1], segA[0], segB[0])
    angle2 = cal_angle(segA[-1], segA[1], segB[0])
    angle3 = cal_angle(segA[-1], segB[-2], segB[0])
    angle = 0.8 * angle1 + 0.1 * angle2 + 0.1 * angle3
    i = 0
    if graph.degree[segA[-1]] == 1:
        i += 1
    angle = (0.8 ** i) * angle
    return [angle, segA[1]]

def segs_angle_othertwosege(segA, segB):
    angle1 = cal_angle(segA[-1], segA[0], segB[-1])
    angle2 = cal_angle(segA[-1], segA[1], segB[-1])
    angle3 = cal_angle(segA[-1], segB[1], segB[-1])
    angle = 0.8 * angle1 + 0.1 * angle2 + 0.1 * angle3
    return [angle, segA[1]]
def segs_angle_othertwosege_new(segA, segB, graph):
    angle1 = cal_angle(segA[-1], segA[0], segB[-1])
    angle2 = cal_angle(segA[-1], segA[1], segB[-1])
    angle3 = cal_angle(segA[-1], segB[1], segB[-1])
    angle = 0.8 * angle1 + 0.1 * angle2 + 0.1 * angle3
    i = 0
    if graph.degree[segA[-1]] == 1:
        i += 1
    if graph.degree[segB[-1]] == 1:
        i += 1
    angle = (0.8 ** i) * angle
    return [angle, segA[1]]
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
            for edge in list(graph.edges(nbunch=b)):
                if edge[1] not in seg:
                    if edge[1] not in Disk_node:
                        seg.append(edge[1])
                    if edge[1] in Disk_node:
                        m = 1
                        if edge[0] not in close_disc:
                            close_disc.append(edge[0])
            if m != 0:
                break

            for edge in list(graph.edges(nbunch=seg[-1])):
                if edge[1] in Disk_node:
                    if edge[0] not in close_disc:
                        close_disc.append(edge[0])
                    break
                break
            if seg[-1] in B_g4 or edge[1] in Disk_node or seg[-1] == b:
                break
        return seg, close_disc

def find_path_new(seg):
    len_seg = len(seg)
    if len_seg < 5:
        path_G_direct = seg
    elif len_seg >= 5 and len_seg <= 7:
        path_G_direct = [seg[0], seg[1], seg[-2], seg[-1]]
    elif len_seg > 7 and len_seg < 13:
        path_G_direct = [seg[0], seg[1], seg[int(len(seg) / 2)], seg[-2], seg[-1]]
    else:
        N_seg = int(len_seg / 6)
        N_interval = round(len_seg / N_seg)
        path_G_direct = []
        path_G_direct.append(seg[0])
        path_G_direct.append(seg[1])
        for l4 in range(1, N_seg + 1):
            if l4 * N_interval >= len_seg - 2:
                break
            path_G_direct.append(seg[l4 * N_interval])
        path_G_direct.append(seg[-2])
        path_G_direct.append(seg[-1])
    return path_G_direct

def trace_source(g4, Disk_node, H, Disc_size, D_t_x, image_label):
    branchnodes_list = {}
    G_direct = nx.DiGraph()
    j = 0
    close_disc = []
    g4_copy = g4.copy()
    nrb_1_g4 = []
    nrb_bn_g4 = []
    for node in g4:
        if g4.degree[node] == 1:
            nrb_1_g4.append(node)
        if g4.degree[node] > 2:
            nrb_bn_g4.append(node)
    B_g4_copy = nrb_1_g4 + nrb_bn_g4

    i_c = 0
    node_correct = []
    for c in sorted(nx.connected_components(g4), key=len, reverse=True):
        i_c += 1
        D_g4 = []
        B_g4 = []
        c = g4.subgraph(c).copy()
        for node in c.nodes:
            if g4.degree[node] == 1 and node not in Disk_node:
                D_g4.append(node)
            if g4.degree[node] > 2 and node not in Disk_node:
                B_g4.append(node)
        D_g4_copy = D_g4.copy()
        if len(D_g4) == 0:
            break
        if len(c.nodes)/len(D_g4) < 15:
            continue
        if i_c > 1:
            if i_c == 2:
                num = len(c.nodes)
            if len(c.nodes)/num < 1/2:
                continue
            j_c = 0
            for node in c.nodes:
                if (image_label[node[0], node[1], :] != [0, 0, 0]).any():
                    j_c += 1
            if j_c/len(c.nodes) < 0.5:
                continue

        '''20220525添加'''
        D_g4_1 = []
        D_g4_2 = []
        for node in D_g4:
            if node in D_t_x:
                D_g4_1.append(node)
            else:
                D_g4_2.append(node)
        D_g4 = D_g4_1 + D_g4_2

        D_g4_1 = []
        while len(D_g4) > 0:
            d = D_g4.pop(-1)
            check_last_seg = [0, 0]
            path_G_direct = [0]
            check_penultimate_seg = []
            nrb = 0
            path_last = 0
            while 1:
                node = d
                if g4.degree[node] == 1:
                    if node in G_direct:
                        break
                    edge = list(g4.edges(nbunch=node))
                    nrb = edge[0][1]
                    seg, close_disc = find_seg_new(node, nrb, close_disc, g4, Disk_node, B_g4_copy)
                    branchnodes_list[f'{node},{nrb}'] = seg
                    if len(seg) == 2:
                        g4.remove_node(node)
                        break
                    path_G_direct = [seg[0], seg[int(len(seg) / 2)], seg[-1]]

                    dist1 = math.sqrt((node[0] - H[0]) ** 2 + (node[1] - H[1]) ** 2)
                    dist2 = math.sqrt((path_G_direct[-1][0] - H[0]) ** 2 + (path_G_direct[-1][1] - H[1]) ** 2)
                    abs_dist = math.sqrt((path_G_direct[-1][0] - node[0]) ** 2 + (path_G_direct[-1][1] - node[1]) ** 2)

                    if (dist2 - dist1) / abs_dist > 0.8 and (dist2 < Disc_size * 8 or dist1 < Disc_size * 4):
                        D_g4_1.append(node)
                        break

                    if dist1 < Disc_size * 3:
                        if (abs_dist < 40 and (dist2 - dist1) / abs_dist > 0.8) or (abs_dist >= 40 and (dist2 - dist1) / abs_dist > 0.6):
                            D_g4_1.append(node)
                            break

                    nx.add_path(G_direct, path_G_direct)
                    path_last = path_G_direct
                    check_last_seg = [node, nrb]
                    d = seg[-1]
                    continue
                if node in close_disc:
                    break
                if g4_copy.degree[node] > 2:
                    dist = math.sqrt((node[0] - H[0]) ** 2 + (node[1] - H[1]) ** 2)
                    if g4.degree[node] == 4 and dist < 3 * Disc_size:
                        break
                    nrb_new = (0, 0)
                    i = 0
                    nrb_list = []
                    [angle, nrb] = [0, 0]
                    [angle_, nrb] = [0, 0]
                    bili = 0
                    degree_node = g4.degree[node]
                    for edge in list(g4.edges(nbunch=node)):
                        nrb = edge[1]
                        seg, close_disc = find_seg_new(node, nrb, close_disc, g4, Disk_node, B_g4_copy)
                        branchnodes_list[f'{node},{nrb}'] = seg
                        if seg[-1] != path_G_direct[0]:
                            if seg[-1] in B_g4 or seg[-1] in close_disc or seg[-1] in D_g4 or seg[-1] in D_g4_1:
                                segA_0 = branchnodes_list[f'{node},{nrb}']
                                segA = [segA_0[0], segA_0[1], segA_0[int(len(segA_0) / 2)], segA_0[-2], segA_0[-1]]
                                segB_0 = branchnodes_list[f'{check_last_seg[0]},{check_last_seg[1]}']
                                segB = [segB_0[0], segB_0[1], segB_0[int(len(segB_0) / 2)], segB_0[-2], segB_0[-1]]
                                [angle0, nbr0] = segs_angle_new(segA, segB, g4)
                                i += 1
                                nrb_list.append(nrb)
                                if angle0 > angle_:
                                    [angle_, nrb_new_] = [angle0, nbr0]
                                dist1 = math.sqrt((segA_0[0][0] - H[0]) ** 2 + (segA_0[0][1] - H[1]) ** 2)
                                dist2 = math.sqrt((segA_0[-1][0] - H[0]) ** 2 + (segA_0[-1][1] - H[1]) ** 2)
                                abs_dist = math.sqrt(
                                    (segA_0[-1][0] - segA_0[0][0]) ** 2 + (segA_0[-1][1] - segA_0[0][1]) ** 2)
                                bili0 = (dist1 - dist2) / (abs_dist + 1e-6)
                                if g4.degree[node] == 3 and g4.degree[segB[0]] == 1:
                                    angle1 = angle0 + bili0 * 90
                                else:
                                    angle1 = angle0 + bili0 * 10
                                if angle1 > angle:
                                    [angle, nrb_new] = [angle1, nbr0]

                    nrb = nrb_new
                    if i == 2:
                        segC_0 = branchnodes_list[f'{node},{nrb_list[0]}']
                        segD_0 = branchnodes_list[f'{node},{nrb_list[1]}']
                        segC = [segC_0[0], segC_0[1], segC_0[int(len(segC_0) / 2)], segC_0[-2], segC_0[-1]]
                        segD = [segD_0[0], segD_0[1], segD_0[int(len(segD_0) / 2)], segD_0[-2], segD_0[-1]]
                        [angle1, nbr1] = segs_angle_othertwosege_new(segC, segD, g4)
                        if angle1 > angle:
                            nrb = nbr1
                            break
                    if angle < 1:
                        break
                    seg_next0 = branchnodes_list[f'{node},{nrb}']
                    seg_next = [seg_next0[0], seg_next0[1], seg_next0[int(len(seg_next0) / 2)], seg_next0[-2], seg_next0[-1]]
                    if i > 1:
                        dist1 = math.sqrt((node[0]-H[0])**2 + (node[1]-H[1])**2)
                        dist2 = math.sqrt((seg_next[-1][0]-H[0])**2 + (seg_next[-1][1]-H[1])**2)
                        abs_dist = math.sqrt((seg_next[-1][0] - node[0]) ** 2 + (seg_next[-1][1] - node[1]) ** 2)
                        if dist1 < Disc_size * 3:
                            if (dist2 - dist1) / abs_dist > 0.5:
                                break
                        abs_dist = math.sqrt((seg_next[-1][0]-node[0])**2 + (seg_next[-1][1]-node[1])**2)
                        if dist2 - dist1 > 0:
                            if abs(seg_next[0][0]-H[0]) > H[0]/3 and abs(seg_next[-1][0]-H[0]) > H[0]/3:
                                if g4.degree[seg_next[0]] == 3 and (dist2 - dist1) / abs_dist > 0.5:
                                    break
                            else:
                                if (dist2 - dist1) / abs_dist > 0.5:
                                    break

                    path_G_direct = [seg_next[0], seg_next[2], seg_next[-1]]

                    dist1 = math.sqrt((seg_next0[0][0] - H[0]) ** 2 + (seg_next0[0][1] - H[1]) ** 2)
                    dist2 = math.sqrt((seg_next0[-1][0] - H[0]) ** 2 + (seg_next0[-1][1] - H[1]) ** 2)
                    abs_dist = math.sqrt(
                        (seg_next0[-1][0] - seg_next0[0][0]) ** 2 + (seg_next0[-1][1] - seg_next0[0][1]) ** 2)
                    bili0 = (dist1 - dist2) / (abs_dist + 1e-6)
                    if bili0 > 0.8 and g4.degree[segB[0]] == 1 and dist1 < Disc_size * 2:
                        if path_G_direct[1] in G_direct:
                            break
                        for node in path_last:
                            G_direct.remove_node(node)
                        break

                    check_penultimate_seg = check_last_seg
                    check_last_seg = [node, nrb]
                    if path_G_direct[1] in G_direct.nodes or seg_next0[-int(len(seg_next0) / 2) - 1] in G_direct.nodes:
                        break
                    dist = math.sqrt((path_G_direct[-1][0] - H[0]) ** 2 + (path_G_direct[-1][1] - H[1]) ** 2)
                    if g4.degree[path_G_direct[-1]] == 1 and dist > Disc_size * 2:
                        break

                    nx.add_path(G_direct, path_G_direct)
                    if g4.degree[seg_next[0]] == 3 and g4.degree[segB[0]] == 1 and nrb_new != nrb_new_:
                        node_correct.append(d)
                    d = seg_next[-1]
                    if d in close_disc:
                        break
                    if g4.degree[d] == 1:
                        break
            j += 1
            if len(D_g4) == 0:
                break
        seg_midpoint_list = {}
        B_g4_copy2 = B_g4_copy.copy()
        node_join = []
        while 1:
            missing1_degree_branchnode = []
            for node in B_g4:
                if node in G_direct.nodes:
                    if node in Disk_node:
                        continue
                    if g4.degree[node] - G_direct.degree[node] == 1:
                        missing1_degree_branchnode.append(node)
            if len(missing1_degree_branchnode) == 0:
                break
            seg_miss = []
            for node in missing1_degree_branchnode:
                if g4.degree[node] - G_direct.degree[node] == 0 and len(seg_miss) != 0:
                    if seg_miss[-1] != node and seg_miss[0] != node:
                        continue

                nrb_miss_inDisk = []
                nrb_miss = []
                seg_in = []
                seg_out = []
                for edge in list(g4.edges(nbunch=node)):
                    nrb_miss.append(edge[1])
                for i in range(0, len(nrb_miss)):
                    if nrb_miss[i] in Disk_node:
                        nrb_miss_inDisk.append(nrb_miss[i])
                        B_g4 = [x for x in B_g4 if x != node]
                        break
                    seg, close_disc = find_seg_new(node, nrb_miss[i], close_disc, g4, Disk_node, B_g4_copy)
                    seg1 = [seg[0], seg[1], seg[int(len(seg) / 2)], seg[-2], seg[-1]]
                    seg2 = [seg[0], seg[1], seg[-int(len(seg) / 2) - 1], seg[-2], seg[-1]]
                    if len(seg) < 3:
                        seg1 = [seg[0], seg[0], seg[-1]]
                        seg2 = [seg[0], seg[0], seg[-1]]
                    if seg1[2] in G_direct.nodes or seg2[2] in G_direct.nodes:
                        if seg1[2] in G_direct.nodes:
                            seg = seg1
                        if seg2[2] in G_direct.nodes:
                            seg = seg2
                        for edge in list(G_direct.in_edges(nbunch=seg[2])):
                            if edge[0] == seg[0]:
                                seg_out.append(seg)
                        for edge in list(G_direct.out_edges(nbunch=seg[2])):
                            if edge[1] == seg[0]:
                                seg_in.append(seg)
                        continue
                    else:
                        seg_miss = seg
                if len(seg_miss) < 3:
                    B_g4 = [x for x in B_g4 if x != node]
                    continue

                if len(nrb_miss_inDisk) != 0:
                    continue

                nextBN_in = []
                nextBN_out = []
                if seg_miss[-1] in G_direct.nodes:
                    if G_direct.degree[seg_miss[-1]] != 0:
                        nextBN_nrb = []
                        for edge in g4.edges(nbunch=seg_miss[-1]):
                            nextBN_nrb.append(edge[1])
                        for l1 in range(0, len(nextBN_nrb)):
                            seg_nextBN, close_disc = find_seg_new(seg_miss[-1], nextBN_nrb[l1], close_disc, g4,
                                                                  Disk_node, B_g4_copy)
                            seg_nextBN_1 = [seg_nextBN[0], seg_nextBN[1], seg_nextBN[int(len(seg_nextBN) / 2)], seg_nextBN[-2], seg_nextBN[-1]]
                            seg_nextBN_2 = [seg_nextBN[0], seg_nextBN[1], seg_nextBN[-int(len(seg_nextBN) / 2) - 1], seg_nextBN[-2], seg_nextBN[-1]]
                            if seg_nextBN_1[2] in G_direct.nodes or seg_nextBN_2[2] in G_direct.nodes:
                                if seg_nextBN_1[2] in G_direct.nodes:
                                    seg_nextBN = seg_nextBN_1
                                if seg_nextBN_2[2] in G_direct.nodes:
                                    seg_nextBN = seg_nextBN_2
                                for edge in list(G_direct.in_edges(nbunch=seg_nextBN[2])):
                                    if edge[0] == seg_miss[-1]:
                                        nextBN_out.append(seg_nextBN)
                                for edge in list(G_direct.out_edges(nbunch=seg_nextBN[2])):
                                    if edge[1] == seg_miss[-1]:
                                        nextBN_in.append(seg_nextBN)
                angle_in = 0
                for seg in seg_in:
                    [angle, nrb_in] = segs_angle_othertwosege(seg, seg_miss)
                    angle_in += angle
                for seg in nextBN_out:
                    [angle, nrb_in] = segs_angle(seg, seg_miss)
                    w = 1 / len(nextBN_nrb)
                    angle_in += w * angle
                angle_out = 0
                for seg in seg_out:
                    [angle, nrb_out] = segs_angle_othertwosege(seg, seg_miss)
                    angle_out += angle
                for seg in nextBN_in:
                    [angle, nrb_out] = segs_angle(seg, seg_miss)
                    w = 1 / len(nextBN_nrb)
                    angle_out += w * angle
                if angle_in > angle_out:
                    path_G_direct = [seg_miss[0], seg_miss[int(len(seg_miss) / 2)], seg_miss[-1]]
                    angle_next = angle_in - angle_out
                else:
                    path_G_direct = [seg_miss[-1], seg_miss[int(len(seg_miss) / 2)], seg_miss[0]]
                    angle_next = angle_out - angle_in
                dist_first = math.sqrt((path_G_direct[0][0] - H[0]) ** 2 + (path_G_direct[0][1] - H[1]) ** 2)
                dist_last = math.sqrt((path_G_direct[-1][0] - H[0]) ** 2 + (path_G_direct[-1][1] - H[1]) ** 2)
                abs_dist = math.sqrt((path_G_direct[-1][0] - path_G_direct[0][0]) ** 2 + (path_G_direct[-1][1] - path_G_direct[0][1]) ** 2)
                if dist_last - dist_first > 0:
                    if dist_first < Disc_size * 3:
                        if (dist_last - dist_first) / abs_dist > abs(angle_in - angle_out) / ((angle_in + angle_out) / 2):
                            path_G_direct = path_G_direct[::-1]
                    else:
                        if (dist_last - dist_first) / abs_dist > 0.5:
                            if (dist_last - dist_first) / abs_dist > abs(angle_in - angle_out) / (
                                    (angle_in + angle_out) / 2):
                                path_G_direct = path_G_direct[::-1]
                            if len(seg_in) == 0:
                                path_G_direct = path_G_direct[::-1]

                score_segmiss_v = 0
                score_segmiss_a = 0
                for node in seg_miss:
                    if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                        score_segmiss_v += 1
                    if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                        score_segmiss_a += 1
                proportion_v_segmiss = score_segmiss_v / (score_segmiss_v + score_segmiss_a + 1e-8)
                proportion_a_segmiss = score_segmiss_a / (score_segmiss_v + score_segmiss_a + 1e-8)
                (i1, j1) = path_G_direct[0]
                score_v1, score_a1 = find_10x10(i1, j1, image_label)
                (i2, j2) = path_G_direct[-1]
                score_v2, score_a2 = find_10x10(i2, j2, image_label)
                dist = math.sqrt((path_G_direct[-1][0] - H[0]) ** 2 + (path_G_direct[-1][1] - H[1]) ** 2)
                nrbout_segmiss0 = []
                for edge in list(G_direct.out_edges(nbunch=seg_miss[0])):
                    nrbout_segmiss0.append(edge[1])
                nrbout_segmiss1 = []
                for edge in list(G_direct.out_edges(nbunch=seg_miss[-1])):
                    nrbout_segmiss1.append(edge[1])
                angle0, i0 = 0, 0
                for node_ in nrbout_segmiss0:
                    seg, close_disc = find_seg_new(seg_miss[0], node_, close_disc, g4, Disk_node, B_g4_copy)
                    (angle, nrb) = segs_angle_othertwosege(seg, seg_miss)
                    angle0 += angle
                    i0 += 1
                angle0 = angle0/(i0 + 1e-8)
                angle1, i1 = 0, 0
                for node_ in nrbout_segmiss1:
                    seg, close_disc = find_seg_new(seg_miss[-1], node_, close_disc, g4, Disk_node, B_g4_copy)
                    (angle, nrb) = segs_angle_othertwosege(seg, seg_miss[::-1])
                    angle1 += angle
                    i1 += 1
                angle1 = angle1 / (i1 + 1e-8)
                if proportion_v_segmiss > proportion_a_segmiss and score_v2 / (score_v2 + score_a2 + 1e-6) < score_v1 / (
                        score_v1 + score_a1 + 1e-6) and score_a2 / (score_v2 + score_a2 + 1e-6) > 0.4 and dist > Disc_size * 2 and score_a1 == 0 and abs(angle1 - angle0) < 40:
                    path_G_direct = path_G_direct[::-1]
                    if path_G_direct[0] not in node_join:
                        node_join.append(path_G_direct[0])
                if proportion_a_segmiss > proportion_v_segmiss and score_a2 / (score_v2 + score_a2 + 1e-6) < score_a1 / (
                        score_v1 + score_a1 + 1e-6) and score_v2 / (
                        score_v2 + score_a2 + 1e-6) > 0.4 and dist > Disc_size * 2 and score_v1 == 0 and abs(angle1 - angle0) <40:
                    path_G_direct = path_G_direct[::-1]
                    if path_G_direct[0] not in node_join:
                        node_join.append(path_G_direct[0])

                if path_G_direct[1] in G_direct.nodes:
                    if list(G_direct.in_edges(nbunch=path_G_direct[1]))[0][0] == path_G_direct[0]:
                        break
                    if list(G_direct.out_edges(nbunch=path_G_direct[1]))[0][1] == path_G_direct[0]:
                        angle_next_before = seg_midpoint_list[f'{path_G_direct[1]}']
                        if angle_next < angle_next_before[0]:
                            break
                        if angle_next > angle_next_before[0]:
                            G_direct.remove_edge(path_G_direct[-1], path_G_direct[1])
                            G_direct.remove_edge(path_G_direct[1], path_G_direct[0])
                if seg_miss[-int(len(seg_miss) / 2) - 1] != path_G_direct[1]:
                    if seg_miss[-int(len(seg_miss) / 2) - 1] in G_direct.nodes:
                        angle_next_before = seg_midpoint_list[f'{path_G_direct[1]}']
                        if angle_next < angle_next_before[0]:
                            break
                        if angle_next > angle_next_before[0]:
                            G_direct.remove_edge(path_G_direct[-1], seg_miss[-int(len(seg_miss) / 2) - 1])
                            G_direct.remove_edge(seg_miss[-int(len(seg_miss) / 2) - 1], path_G_direct[0])
                nx.add_path(G_direct, path_G_direct)
                seg_midpoint_list[f'{path_G_direct[1]}'] = [angle_next]

                d = path_G_direct[-1]
                if seg_miss[0] == path_G_direct[0]:
                    seg_last = seg_miss
                if seg_miss[0] == path_G_direct[-1]:
                    seg_last = seg_miss[::-1]
                while 1:
                    node = d
                    if g4.degree[node] - G_direct.degree[node] > 0:
                        [angle, nrb] = [0, 0]
                        i_segmiss = 0
                        for edge in list(g4.edges(nbunch=node)):
                            nrb = edge[1]
                            seg, close_disc = find_seg_new(node, nrb, close_disc, g4, Disk_node, B_g4_copy)
                            branchnodes_list[f'{node},{nrb}'] = seg
                            if seg[-1] != seg_last[0]:
                                seg_0 = [seg[0], seg[1], seg[int(len(seg) / 2)], seg[-2], seg[-1]]
                                [angle0, nbr0] = segs_angle_new(seg_0, seg_last, g4)
                                i_segmiss += 1
                                if angle0 > angle:
                                    [angle, nrb_new] = [angle0, nbr0]
                        nrb = nrb_new
                        seg_next0 = branchnodes_list[f'{node},{nrb}']
                        seg_next = [seg_next0[0], seg_next0[1], seg_next0[int(len(seg_next0) / 2)], seg_next0[-2],
                                    seg_next0[-1]]

                        if i_segmiss > 1:
                            dist1 = math.sqrt((node[0] - H[0]) ** 2 + (node[1] - H[1]) ** 2)
                            dist2 = math.sqrt((seg_next[-1][0] - H[0]) ** 2 + (seg_next[-1][1] - H[1]) ** 2)
                            abs_dist = math.sqrt((seg_next[-1][0] - node[0]) ** 2 + (seg_next[-1][1] - node[1]) ** 2)
                            if dist2 - dist1 > 0:
                                if (dist2 - dist1) / abs_dist > 0.5:
                                    break

                            in_degree = (0 if len({G_direct.in_degree(nbunch=seg_next[-1])}) == 0 else G_direct.in_degree(nbunch=seg_next[-1]))
                            out_degree = (0 if len({G_direct.out_degree(nbunch=seg_next[-1])}) == 0 else G_direct.out_degree(nbunch=seg_next[-1]))
                            nrb = []
                            for edge in list(g4.edges(nbunch=seg_next[-1])):
                                if edge[1] != seg_next[-2]:
                                    nrb.append(edge[1])
                            Seg = []
                            for node in nrb:
                                seg = find_seg_new(seg_next[-1], node, close_disc, g4, Disk_node, B_g4_copy)
                                Seg.append(seg)
                            seg_in = []
                            for edge in list(G_direct.in_edges(nbunch=seg_next[-1])):
                                for i in range(len(Seg)):
                                    if edge[0] in Seg[i]:
                                        seg_in.append(Seg[i])
                            seg_out = []
                            for edge in list(G_direct.out_edges(nbunch=seg_next[-1])):
                                for i in range(len(Seg)):
                                    if edge[1] in Seg[i]:
                                        seg_out.append(Seg[i])
                            if in_degree != 0:
                                (angle, nrb) = (180, 0)
                                angle_in_min = 0
                                for seg in seg_in:
                                    (angle0, nrb) = segs_angle(seg_next, seg)
                                    if angle0 < angle:
                                        angle = angle0
                                    angle_in_min = angle
                            if out_degree != 0:
                                (angle, nrb) = (0, 0)
                                angle_out_max = 180
                                for seg in seg_out:
                                    (angle0, nrb) = segs_angle(seg_next, seg)
                                    if angle0 > angle:
                                        angle = angle0
                                    angle_out_max = angle
                            if (in_degree != 0 and angle_in_min > 110) and (out_degree != 0 and angle_out_max < 90):
                                break
                            if in_degree != 0 and angle_in_min > 110 and out_degree == 0:
                                break
                            if out_degree != 0 and angle_out_max < 90 and in_degree == 0:
                                break

                            if abs(seg_next[0][0] - H[0]) > H[0] / 3 and abs(seg_next[-1][0] - H[0]) > H[0] / 3:
                                abs_dist = math.sqrt(
                                    (seg_next[-1][0] - node[0]) ** 2 + (seg_next[-1][1] - node[1]) ** 2)
                                if dist2 - dist1 > 0 and (dist2 - dist1) / abs_dist > 0.5:
                                    break
                            if abs(seg_next[0][0] - H[0]) < H[0] / 5 and abs(seg_next[-1][0] - H[0]) < H[0] / 5:
                                abs_dist = math.sqrt(
                                    (seg_next[-1][0] - node[0]) ** 2 + (seg_next[-1][1] - node[1]) ** 2)
                                if dist2 - dist1 > 0 and (dist2 - dist1) / abs_dist > 0.5:
                                    break

                        path_G_direct = [seg_next[0], seg_next[2], seg_next[-1]]
                        if path_G_direct[1] in G_direct.nodes or seg_next0[-int(len(seg_next0) / 2) - 1] in G_direct.nodes:
                            break
                        nx.add_path(G_direct, path_G_direct)
                        seg_last = seg_next
                        d = seg_next[-1]
                        if d in close_disc:
                            break
                    else:
                        break

        missing1_degree_branchnode1 = []
        for node in B_g4_copy2:
            if node in G_direct.nodes:
                if node in Disk_node:
                    continue
                if g4.degree[node] - G_direct.degree[node] > 0:
                    missing1_degree_branchnode1.append(node)
        for node in missing1_degree_branchnode1:
            if g4.degree[node] - G_direct.degree[node] == 0 and len(seg_miss) != 0:
                if seg_miss[-1] != node and seg_miss[0] != node:
                    continue
            nrb_miss = []
            for edge in list(g4.edges(nbunch=node)):
                nrb_miss.append(edge[1])
            for i in range(0, len(nrb_miss)):
                if nrb_miss[i] in Disk_node:
                    continue
                seg, close_disc = find_seg_new(node, nrb_miss[i], close_disc, g4, Disk_node, B_g4_copy)
                dist_first = math.sqrt((seg[0][0] - H[0]) ** 2 + (seg[0][1] - H[1]) ** 2)
                dist_last = math.sqrt((seg[-1][0] - H[0]) ** 2 + (seg[-1][1] - H[1]) ** 2)
                if len(seg) < 8:
                    continue
                elif seg[int(len(seg) / 2)] in G_direct.nodes:
                    continue
                elif seg[-int(len(seg) / 2) - 1] in G_direct.nodes:
                    continue
                else:
                    if dist_first > dist_last:
                        path_G_direct = [seg[0], seg[int(len(seg) / 2)], seg[-1]]
                    else:
                        path_G_direct = [seg[-1], seg[int(len(seg) / 2)], seg[0]]
                nx.add_path(G_direct, path_G_direct)

        Bn_false_3 = []
        Bn_false_4_1 = []
        Bn_false_4_2 = []
        for node in B_g4:
            in_degree = (0 if len({G_direct.in_degree(nbunch=node)}) == 0 else G_direct.in_degree(nbunch=node))
            out_degree = (0 if len({G_direct.out_degree(nbunch=node)}) == 0 else G_direct.out_degree(nbunch=node))
            if g4.degree[node] == 3:
                if in_degree == 0 or out_degree == 0:
                    Bn_false_3.append(node)
            if g4.degree[node] == 4:
                if in_degree == 0 or out_degree == 0:
                    Bn_false_4_1.append(node)
                if in_degree == 1 or out_degree == 1:
                    Bn_false_4_2.append(node)

        for node in Bn_false_3:
            Bn_false_3_list = {}
            Bn_false_3_nrb = []
            angle_twoseg = {}
            for edge in list(g4.edges(nbunch=node)):
                seg, close_disc = find_seg_new(node, edge[1], close_disc, g4, Disk_node, B_g4_copy)
                Bn_false_3_list[f'{node},{edge[1]}'] = seg
                Bn_false_3_nrb.append(edge[1])
            angle2 = 360
            for l3 in range(0, len(Bn_false_3_nrb)-1):
                [angle, nrb_in] = segs_angle_othertwosege(Bn_false_3_list[f'{node},{Bn_false_3_nrb[l3]}'], Bn_false_3_list[f'{node},{Bn_false_3_nrb[l3+1]}'])
                angle_twoseg[f'{angle}'] = [Bn_false_3_nrb[l3], Bn_false_3_nrb[l3+1]]
                if angle < angle2:
                    angle2 = angle
            for node1 in Bn_false_3_nrb:
                if node1 not in angle_twoseg[f'{angle2}']:
                    seg_false = Bn_false_3_list[f'{node},{node1}']
            if G_direct.degree(nbunch=node) < 3:
                continue
            if seg_false[-int(len(seg_false) / 2) - 1] in G_direct:
                seg_false_midpoint = seg_false[-int(len(seg_false) / 2) - 1]
            if seg_false[int(len(seg_false) / 2)] in G_direct:
                seg_false_midpoint = seg_false[int(len(seg_false) / 2)]
            if list(G_direct.in_edges(nbunch=seg_false_midpoint))[0][0] == seg_false[-1]:
                G_direct.remove_edge(seg_false[-1], seg_false_midpoint)
                G_direct.remove_edge(seg_false_midpoint, seg_false[0])
                nx.add_path(G_direct, [seg_false[0], seg_false_midpoint, seg_false[-1]])
                continue
            if list(G_direct.out_edges(nbunch=seg_false_midpoint))[0][1] == seg_false[-1]:
                G_direct.remove_edge(seg_false[0], seg_false_midpoint)
                G_direct.remove_edge(seg_false_midpoint, seg_false[-1])
                nx.add_path(G_direct, [seg_false[-1], seg_false_midpoint, seg_false[0]])

        for node in Bn_false_4_2:
            dist1 = math.sqrt((node[0] - H[0]) ** 2 + (node[1] - H[1]) ** 2)
            if dist1 < Disc_size * 3:
                continue


    Bn_3_false = []
    B_g4_copy1 = B_g4_copy.copy()
    G_direct_old = G_direct.copy()
    seg_false_list = []
    for node in B_g4_copy1:
        if node in node_join:
            B_g4_copy1 = [x for x in B_g4_copy1 if x != node]
            continue
        dist = math.sqrt((node[0] - H[0]) ** 2 + (node[1] - H[1]) ** 2)
        if dist < Disc_size * 2:
            B_g4_copy1 = [x for x in B_g4_copy1 if x != node]
            continue
        if node in node_correct:
            a, v = 0, 0
            for edge in list(g4.edges(nbunch=node)):
                seg, close_disc = find_seg_new(branchnode=node, nbr=edge[1], close_disc=close_disc, graph=g4,
                                                                    Disk_node=Disk_node, B_g4=B_g4_copy)
                score_segmiss_v = 0
                score_segmiss_a = 0
                for node_ in seg:
                    if (image_label[node_[0], node_[1], :] == [255, 0, 0]).all():
                        score_segmiss_v += 1
                    if (image_label[node_[0], node_[1], :] == [0, 0, 255]).all():
                        score_segmiss_a += 1
                proportion_v_segmiss = score_segmiss_v / (score_segmiss_v + score_segmiss_a + 1e-8)
                proportion_a_segmiss = score_segmiss_a / (score_segmiss_v + score_segmiss_a + 1e-8)
                if proportion_v_segmiss > proportion_a_segmiss:
                    v += 1
                if proportion_v_segmiss < proportion_a_segmiss:
                    a += 1
            if a == 3 or v == 3:
                B_g4_copy1 = [x for x in B_g4_copy1 if x != node]
                continue
        if node not in G_direct:
            continue
        if G_direct.degree[node] == 3:
            nrb_in_B3_ = []
            nrb_out_B3_ = []
            for edge in list(G_direct.in_edges(nbunch=node)):
                nrb_in_B3_.append(edge[0])
            for edge in list(G_direct.out_edges(nbunch=node)):
                nrb_out_B3_.append(edge[1])
            seg_in = []
            seg_out = []
            for edge in list(g4.edges(nbunch=node)):
                seg, close_disc = find_seg_new(branchnode=node, nbr=edge[1], close_disc=close_disc, graph=g4,
                                                                    Disk_node=Disk_node, B_g4=B_g4_copy)
                if seg[-int(len(seg) / 2) - 1] in nrb_in_B3_ or seg[int(len(seg) / 2)] in nrb_in_B3_:
                    seg_in.append(seg)
                if seg[-int(len(seg) / 2) - 1] in nrb_out_B3_ or seg[int(len(seg) / 2)] in nrb_out_B3_:
                    seg_out.append(seg)
            seg_false = []
            if len(nrb_in_B3_) == 2:
                seg_in_1 = seg_in[0]
                seg_in_1_ = [seg_in_1[0], seg_in_1[1], seg_in_1[int(len(seg_in_1) / 2)], seg_in_1[-2], seg_in_1[-1]]
                seg_in_2 = seg_in[1]
                seg_in_2_ = [seg_in_2[0], seg_in_2[1], seg_in_2[int(len(seg_in_2) / 2)], seg_in_2[-2], seg_in_2[-1]]
                (angle, nrb) = segs_angle_othertwosege(seg_in_1_, seg_in_2_)
                if angle > 135:
                    Bn_3_false.append(node)
                    B_g4_copy1.append(node)
                    if G_direct.degree[seg_in_1[-1]] == 1:
                        if G_direct.degree[seg_in_2[-1]] == 1:
                            B_g4_copy1 = [x for x in B_g4_copy1 if x != node]
                            continue
                        seg_out_ = [seg_out[0][0], seg_out[0][1], seg_out[0][int(len(seg_out[0]) / 2)], seg_out[0][-2],
                                    seg_out[0][-1]]
                        (angle1, nrb1) = segs_angle_othertwosege(seg_out_, seg_in_2_)
                        if angle1 > 130:
                            node_outlast = seg_out[0][-1]
                            node_in2last = seg_in_2[-1]
                            nrbout_outlast = []
                            nrbin_outlast = []
                            segin_outlast = []
                            segout_outlast = []
                            for edge in list(G_direct.in_edges(nbunch=node_outlast)):
                                if edge[0] not in seg_out[0]:
                                    nrbin_outlast.append(edge[0])
                            for edge in list(G_direct.out_edges(nbunch=node_outlast)):
                                nrbout_outlast.append(edge[1])
                            for edge in list(g4.edges(nbunch=node_outlast)):
                                seg, close_disc = find_seg_new(branchnode=node_outlast, nbr=edge[1],
                                                               close_disc=close_disc,
                                                               graph=g4,
                                                               Disk_node=Disk_node, B_g4=B_g4_copy)
                                if seg[-int(len(seg) / 2) - 1] in nrbin_outlast or seg[
                                    int(len(seg) / 2)] in nrbin_outlast:
                                    segin_outlast.append(seg)
                                if seg[-int(len(seg) / 2) - 1] in nrbout_outlast or seg[
                                    int(len(seg) / 2)] in nrbout_outlast:
                                    segout_outlast.append(seg)
                            angle3 = 0
                            for seg in segin_outlast:
                                seg_ = [seg[0], seg[1], seg[int(len(seg) / 2)], seg[-2], seg[-1]]
                                (angle1, nrb1) = segs_angle_othertwosege(seg_out_[::-1], seg_)
                                angle3 += (90 - angle1)
                            for seg in segout_outlast:
                                seg_ = [seg[0], seg[1], seg[int(len(seg) / 2)], seg[-2], seg[-1]]
                                (angle1, nrb1) = segs_angle_othertwosege(seg_out_[::-1], seg_)
                                angle3 += (angle1 - 90)
                            nrbout_in2last = []
                            nrbin_in2last = []
                            segin_in2last = []
                            segout_in2last = []
                            for edge in list(G_direct.in_edges(nbunch=node_in2last)):
                                nrbin_in2last.append(edge[0])
                            for edge in list(G_direct.out_edges(nbunch=node_in2last)):
                                if edge[1] not in seg_in_2:
                                    nrbout_in2last.append(edge[1])
                            for edge in list(g4.edges(nbunch=node_in2last)):
                                seg, close_disc = find_seg_new(branchnode=node_in2last, nbr=edge[1],
                                                               close_disc=close_disc,
                                                               graph=g4,
                                                               Disk_node=Disk_node, B_g4=B_g4_copy)
                                if seg[-int(len(seg) / 2) - 1] in nrbin_in2last or seg[
                                    int(len(seg) / 2)] in nrbin_in2last:
                                    segin_in2last.append(seg)
                                if seg[-int(len(seg) / 2) - 1] in nrbout_in2last or seg[
                                    int(len(seg) / 2)] in nrbout_in2last:
                                    segout_in2last.append(seg)
                            angle4 = 0
                            for seg in segin_in2last:
                                seg_ = [seg[0], seg[1], seg[int(len(seg) / 2)], seg[-2], seg[-1]]
                                (angle1, nrb1) = segs_angle_othertwosege(seg_in_2_[::-1], seg_)
                                angle4 += (angle1 - 90)
                            for seg in segout_in2last:
                                seg_ = [seg[0], seg[1], seg[int(len(seg) / 2)], seg[-2], seg[-1]]
                                (angle1, nrb1) = segs_angle_othertwosege(seg_in_2_[::-1], seg_)
                                angle4 += (90 - angle1)
                            if angle3 > 0 and angle4 > 0:
                                B_g4_copy1 = [x for x in B_g4_copy1 if x != node]
                                continue
                        seg_false = seg_in_2
                        if seg_false in seg_false_list:
                            seg_false = seg_in_1
                        seg_false_list.append(seg_false)
                    elif G_direct.degree[seg_in_2[-1]] == 1:
                        seg_out_ = [seg_out[0][0], seg_out[0][1], seg_out[0][int(len(seg_out[0]) / 2)], seg_out[0][-2],
                                    seg_out[0][-1]]
                        (angle1, nrb1) = segs_angle_othertwosege(seg_out_, seg_in_1_)
                        if angle1 > 130:
                            node_outlast = seg_out[0][-1]
                            node_in1last = seg_in_1[-1]
                            nrbout_outlast = []
                            nrbin_outlast = []
                            segin_outlast = []
                            segout_outlast = []
                            for edge in list(G_direct.in_edges(nbunch=node_outlast)):
                                if edge[0] != seg_out[0][-2]:
                                    nrbin_outlast.append(edge[0])
                            for edge in list(G_direct.out_edges(nbunch=node_outlast)):
                                nrbout_outlast.append(edge[1])
                            for edge in list(g4.edges(nbunch=node_outlast)):
                                seg, close_disc = find_seg_new(branchnode=node_outlast, nbr=edge[1],
                                                               close_disc=close_disc,
                                                               graph=g4,
                                                               Disk_node=Disk_node, B_g4=B_g4_copy)
                                if seg[-int(len(seg) / 2) - 1] in nrbin_outlast or seg[
                                    int(len(seg) / 2)] in nrbin_outlast:
                                    segin_outlast.append(seg)
                                if seg[-int(len(seg) / 2) - 1] in nrbout_outlast or seg[
                                    int(len(seg) / 2)] in nrbout_outlast:
                                    segout_outlast.append(seg)
                            angle3 = 0
                            for seg in segin_outlast:
                                seg_ = [seg[0], seg[1], seg[int(len(seg) / 2)], seg[-2], seg[-1]]
                                (angle1, nrb1) = segs_angle_othertwosege(seg_out_[::-1], seg_)
                                angle3 += (90 - angle1)
                            for seg in segout_outlast:
                                seg_ = [seg[0], seg[1], seg[int(len(seg) / 2)], seg[-2], seg[-1]]
                                (angle1, nrb1) = segs_angle_othertwosege(seg_out_[::-1], seg_)
                                angle3 += (angle1 - 90)
                            nrbout_in1last = []
                            nrbin_in1last = []
                            segin_in1last = []
                            segout_in1last = []
                            for edge in list(G_direct.in_edges(nbunch=node_in1last)):
                                nrbin_in1last.append(edge[0])
                            for edge in list(G_direct.out_edges(nbunch=node_in1last)):
                                if edge[1] != seg_in_1[-2]:
                                    nrbout_in1last.append(edge[1])
                            for edge in list(g4.edges(nbunch=node_in1last)):
                                seg, close_disc = find_seg_new(branchnode=node_in1last, nbr=edge[1],
                                                               close_disc=close_disc,
                                                               graph=g4,
                                                               Disk_node=Disk_node, B_g4=B_g4_copy)
                                if seg[-int(len(seg) / 2) - 1] in nrbin_in1last or seg[
                                    int(len(seg) / 2)] in nrbin_in1last:
                                    segin_in1last.append(seg)
                                if seg[-int(len(seg) / 2) - 1] in nrbout_in1last or seg[
                                    int(len(seg) / 2)] in nrbout_in1last:
                                    segout_in1last.append(seg)
                            angle4 = 0
                            for seg in segin_in1last:
                                seg_ = [seg[0], seg[1], seg[int(len(seg) / 2)], seg[-2], seg[-1]]
                                (angle1, nrb1) = segs_angle_othertwosege(seg_in_1_[::-1], seg_)
                                angle4 += (90 - angle1)
                            for seg in segout_in1last:
                                seg_ = [seg[0], seg[1], seg[int(len(seg) / 2)], seg[-2], seg[-1]]
                                (angle1, nrb1) = segs_angle_othertwosege(seg_in_1_[::-1], seg_)
                                angle4 += (angle1 - 90)
                            if angle3 > 0 and angle4 > 0:
                                B_g4_copy1 = [x for x in B_g4_copy1 if x != node]
                                continue
                        seg_false = seg_in_2
                        if seg_false in seg_false_list:
                            seg_false = seg_in_1
                        seg_false_list.append(seg_false)
                    else:
                        seg_in1_lastpoint = seg_in_1[-1]
                        nrb_in1_lastpoint_ = []
                        for edge in list(G_direct.in_edges(nbunch=seg_in1_lastpoint)):
                            nrb_in1_lastpoint_.append(edge[0])
                        seg_in1lastpoint_in = []
                        for edge in list(g4.edges(nbunch=seg_in1_lastpoint)):
                            seg, close_disc = find_seg_new(branchnode=seg_in1_lastpoint, nbr=edge[1], close_disc=close_disc, graph=g4,
                                               Disk_node=Disk_node, B_g4=B_g4_copy)
                            if seg[-int(len(seg) / 2) - 1] in nrb_in1_lastpoint_ or seg[int(len(seg) / 2)] in nrb_in1_lastpoint_:
                                seg_in1lastpoint_in.append(seg)
                        angle_in1 = 180
                        for i in range(0, len(seg_in1lastpoint_in)):
                            seg = seg_in1lastpoint_in[i]
                            (angle1, nrb) = segs_angle(seg, seg_in_1)
                            if angle1 < angle_in1:
                                angle_in1 = angle1
                        seg_in2_lastpoint = seg_in_2[-1]
                        nrb_in2_lastpoint_ = []
                        for edge in list(G_direct.in_edges(nbunch=seg_in2_lastpoint)):
                            nrb_in2_lastpoint_.append(edge[0])
                        seg_in2lastpoint_in = []
                        for edge in list(g4.edges(nbunch=seg_in2_lastpoint)):
                            seg, close_disc = find_seg_new(branchnode=seg_in2_lastpoint, nbr=edge[1], close_disc=close_disc, graph=g4,
                                               Disk_node=Disk_node, B_g4=B_g4_copy)
                            if seg[-int(len(seg) / 2) - 1] in nrb_in2_lastpoint_ or seg[
                                int(len(seg) / 2)] in nrb_in2_lastpoint_:
                                seg_in2lastpoint_in.append(seg)
                        angle_in2 = 180
                        for i in range(0, len(seg_in2lastpoint_in)):
                            seg = seg_in2lastpoint_in[i]
                            (angle1, nrb) = segs_angle(seg, seg_in_2)
                            if angle1 < angle_in2:
                                angle_in2 = angle1
                        if angle_in1 < angle_in2:
                            seg_false = seg_in_1
                            if seg_false in seg_false_list:
                                seg_false = seg_in_2
                        else:
                            seg_false = seg_in_2
                            if seg_false in seg_false_list:
                                seg_false = seg_in_1
                        seg_false_list.append(seg_false)
                if len(seg_false) != 0:
                    if seg_false[-int(len(seg_false) / 2) - 1] in G_direct:
                        seg_false_midpoint = seg_false[-int(len(seg_false) / 2) - 1]
                    if seg_false[int(len(seg_false) / 2)] in G_direct:
                        seg_false_midpoint = seg_false[int(len(seg_false) / 2)]
                    G_direct.remove_edge(seg_false[-1], seg_false_midpoint)
                    G_direct.remove_edge(seg_false_midpoint, seg_false[0])
                    nx.add_path(G_direct, [seg_false[0], seg_false_midpoint, seg_false[-1]])
                    mid_point_in = []
                    mid_point_out = []
                    seg_false_in = []
                    seg_false_out = []
                    for edge in G_direct.in_edges(nbunch=seg_false[-1]):
                        mid_point_in.append(edge[0])
                    for edge in G_direct.out_edges(nbunch=seg_false[-1]):
                        mid_point_out.append(edge[1])
                    for edge in g4.edges(nbunch=seg_false[-1]):
                        seg, close_disc = find_seg_new(branchnode=seg_false[-1], nbr=edge[1], close_disc=close_disc,
                                                       graph=g4, Disk_node=Disk_node, B_g4=B_g4_copy)
                        if seg[::-1] == seg_false:
                            continue
                        if seg[-int(len(seg) / 2) - 1] in mid_point_in or seg[
                            int(len(seg) / 2)] in mid_point_in:
                            seg_false_in.append(seg)
                        if seg[-int(len(seg) / 2) - 1] in mid_point_out or seg[
                            int(len(seg) / 2)] in mid_point_out:
                            seg_false_out.append(seg)
                    angle_falsein = 0
                    for i in range(0, len(seg_false_in)):
                        seg_last = seg_false_in[i]
                        (angle1, nrb) = segs_angle(seg_last, seg_false)
                        if angle1 > angle_falsein:
                            angle_falsein = angle1
                            seg_last_false0 = seg_last
                    seg_last_false = []
                    if angle_falsein > 120:
                        seg_last_false = seg_last_false0
                    if len(seg_last_false) != 0 and G_direct.degree[seg_last_false[-1]] > 1:
                        if seg_last_false[-int(len(seg_last_false) / 2) - 1] in G_direct:
                            seg_false_midpoint = seg_last_false[-int(len(seg_last_false) / 2) - 1]
                        if seg_last_false[int(len(seg_last_false) / 2)] in G_direct:
                            seg_false_midpoint = seg_last_false[int(len(seg_last_false) / 2)]
                        G_direct.remove_edge(seg_last_false[-1], seg_false_midpoint)
                        G_direct.remove_edge(seg_false_midpoint, seg_last_false[0])
                        nx.add_path(G_direct, [seg_last_false[0], seg_false_midpoint, seg_last_false[-1]])
            seg_false = []
            if len(nrb_out_B3_) == 2:
                seg_out_1 = seg_out[0]
                seg_out_1_ = [seg_out_1[0], seg_out_1[1], seg_out_1[int(len(seg_out_1) / 2)], seg_out_1[-2], seg_out_1[-1]]
                seg_out_2 = seg_out[1]
                seg_out_2_ = [seg_out_2[0], seg_out_2[1], seg_out_2[int(len(seg_out_2) / 2)], seg_out_2[-2], seg_out_2[-1]]
                (angle, nrb) = segs_angle_othertwosege(seg_out_1_, seg_out_2_)
                if angle > 135:
                    Bn_3_false.append(node)
                    B_g4_copy1.append(node)
                    if G_direct.degree[seg_out_1[-1]] == 1:
                        seg_false = seg_out_2
                        if seg_false in seg_false_list:
                            seg_false = seg_out_1
                        seg_false_list.append(seg_false)
                    elif G_direct.degree[seg_out_2[-1]] == 1:
                        seg_false = seg_out_1
                        if seg_false in seg_false_list:
                            seg_false = seg_out_2
                        seg_false_list.append(seg_false)
                    else:
                        seg_out1_lastpoint = seg_out_1[-1]
                        nrb_out1_lastpoint_ = []
                        for edge in list(G_direct.out_edges(nbunch=seg_out1_lastpoint)):
                            nrb_out1_lastpoint_.append(edge[1])
                        seg_out1lastpoint_in = []
                        for edge in list(g4.edges(nbunch=seg_out1_lastpoint)):
                            seg, close_disc = find_seg_new(branchnode=seg_out1_lastpoint, nbr=edge[1], close_disc=close_disc, graph=g4,
                                               Disk_node=Disk_node, B_g4=B_g4_copy)
                            if seg[-int(len(seg) / 2) - 1] in nrb_out1_lastpoint_ or seg[
                                int(len(seg) / 2)] in nrb_out1_lastpoint_:
                                seg_out1lastpoint_in.append(seg)
                        angle_out1 = 180
                        for i in range(0, len(seg_out1lastpoint_in)):
                            seg = seg_out1lastpoint_in[i]
                            seg_last1 = []
                            if len(seg) < 5:
                                nrb_last_out_G = []
                                for edge in list(G_direct.out_edges(nbunch=seg[-1])):
                                    nrb_last_out_G.append(edge[1])
                                nrb_last_out1_g = []
                                for edge in list(g4.edges(nbunch=seg[-1])):
                                    nrb_last_out1_g.append(edge[1])
                                seg_last_out1 = []
                                for node_1 in nrb_last_out1_g:
                                    seg_1, close_disc = find_seg_new(branchnode=seg[-1], nbr=node_1,
                                                                     close_disc=close_disc, graph=g4,
                                                                     Disk_node=Disk_node, B_g4=B_g4_copy)
                                    if seg_1[-int(len(seg_1) / 2) - 1] in nrb_last_out_G or seg_1[
                                        int(len(seg_1) / 2)] in nrb_last_out_G:
                                        seg_last_out1.append(seg_1)
                                (angle_1, nrb_1) = (0, 0)
                                for seg_1_ in seg_last_out1:
                                    (angle, nrb) = segs_angle(seg_1, seg)
                                    if angle > angle_1:
                                        (angle_1, nrb_new) = (angle, nrb)
                                seg_last1, close_disc = find_seg_new(branchnode=seg[-1], nbr=nrb_new,
                                                                     close_disc=close_disc, graph=g4,
                                                                     Disk_node=Disk_node, B_g4=B_g4_copy)
                            seg = seg + seg_last1
                            (angle1, nrb) = segs_angle(seg, seg_out_1)
                            if angle1 < angle_out1:
                                angle_out1 = angle1
                        seg_out2_lastpoint = seg_out_2[-1]
                        nrb_out2_lastpoint_ = []
                        for edge in list(G_direct.out_edges(nbunch=seg_out2_lastpoint)):
                            nrb_out2_lastpoint_.append(edge[1])
                        seg_out2lastpoint_in = []
                        for edge in list(g4.edges(nbunch=seg_out2_lastpoint)):
                            seg, close_disc = find_seg_new(branchnode=seg_out2_lastpoint, nbr=edge[1], close_disc=close_disc, graph=g4,
                                               Disk_node=Disk_node, B_g4=B_g4_copy)
                            if seg[-int(len(seg) / 2) - 1] in nrb_out2_lastpoint_ or seg[int(len(seg) / 2)] in nrb_out2_lastpoint_:
                                seg_out2lastpoint_in.append(seg)
                        angle_out2 = 180
                        for i in range(0, len(seg_out2lastpoint_in)):
                            seg = seg_out2lastpoint_in[i]
                            seg_last1 = []
                            if len(seg) < 5:
                                nrb_last_out1_G = []
                                for edge in list(G_direct.out_edges(nbunch=seg[-1])):
                                    nrb_last_out1_G.append(edge[1])
                                nrb_last_out1_g = []
                                for edge in list(g4.edges(nbunch=seg[-1])):
                                    nrb_last_out1_g.append(edge[1])
                                seg_last_out1 = []
                                for node_1 in nrb_last_out1_g:
                                    seg_1, close_disc = find_seg_new(branchnode=seg[-1], nbr=node_1,
                                                                     close_disc=close_disc, graph=g4,
                                                                     Disk_node=Disk_node, B_g4=B_g4_copy)
                                    if seg_1[-int(len(seg_1) / 2) - 1] in nrb_last_out1_G or seg_1[int(len(seg_1) / 2)] in nrb_last_out1_G:
                                        seg_last_out1.append(seg_1)
                                (angle_1, nrb_1) = (0, 0)
                                for seg_1_ in seg_last_out1:
                                    (angle, nrb) = segs_angle(seg_1_, seg)
                                    if angle > angle_1:
                                        (angle_1, nrb_new) = (angle, nrb)
                                seg_last1, close_disc = find_seg_new(branchnode=seg[-1], nbr=nrb_new,
                                                                   close_disc=close_disc, graph=g4,
                                                                   Disk_node=Disk_node, B_g4=B_g4_copy)
                            seg = seg + seg_last1
                            (angle1, nrb) = segs_angle(seg, seg_out_2)
                            if angle1 < angle_out2:
                                angle_out2 = angle1
                        if angle_out1 < angle_out2:
                            seg_false = seg_out_1
                            if seg_false in seg_false_list:
                                seg_false = seg_out_2
                        else:
                            seg_false = seg_out_2
                            if seg_false in seg_false_list:
                                seg_false = seg_out_1
                        seg_false_list.append(seg_false)
                if len(seg_false) != 0:
                    if seg_false[-int(len(seg_false) / 2) - 1] in G_direct:
                        seg_false_midpoint = seg_false[-int(len(seg_false) / 2) - 1]
                    if seg_false[int(len(seg_false) / 2)] in G_direct:
                        seg_false_midpoint = seg_false[int(len(seg_false) / 2)]
                    G_direct.remove_edge(seg_false[0], seg_false_midpoint)
                    G_direct.remove_edge(seg_false_midpoint, seg_false[-1])
                    nx.add_path(G_direct, [seg_false[-1], seg_false_midpoint, seg_false[0]])
                    mid_point_in = []
                    mid_point_out = []
                    seg_false_in = []
                    seg_false_out = []
                    for edge in G_direct.in_edges(nbunch=seg_false[-1]):
                        mid_point_in.append(edge[0])
                    for edge in G_direct.out_edges(nbunch=seg_false[-1]):
                        mid_point_out.append(edge[1])
                    for edge in g4.edges(nbunch=seg_false[-1]):
                        seg, close_disc = find_seg_new(branchnode=seg_false[-1], nbr=edge[1], close_disc=close_disc,
                                                       graph=g4, Disk_node=Disk_node, B_g4=B_g4_copy)
                        if seg[::-1] == seg_false:
                            continue
                        if seg[-int(len(seg) / 2) - 1] in mid_point_in or seg[
                            int(len(seg) / 2)] in mid_point_in:
                            seg_false_in.append(seg)
                        if seg[-int(len(seg) / 2) - 1] in mid_point_out or seg[
                            int(len(seg) / 2)] in mid_point_out:
                            seg_false_out.append(seg)
                    angle_falseout = 0
                    for i in range(0, len(seg_false_out)):
                        seg_last = seg_false_out[i]
                        (angle1, nrb) = segs_angle(seg_last, seg_false)
                        if angle1 > angle_falseout:
                            angle_falseout = angle1
                            seg_last_false0 = seg_last
                    seg_last_false = []
                    if angle_falseout > 120:
                        seg_last_false = seg_last_false0
                    if len(seg_last_false) != 0:
                        if seg_last_false[-int(len(seg_last_false) / 2) - 1] in G_direct:
                            seg_false_midpoint = seg_last_false[-int(len(seg_last_false) / 2) - 1]
                        if seg_last_false[int(len(seg_last_false) / 2)] in G_direct:
                            seg_false_midpoint = seg_last_false[int(len(seg_last_false) / 2)]
                        G_direct.remove_edge(seg_last_false[0], seg_false_midpoint)
                        G_direct.remove_edge(seg_false_midpoint, seg_last_false[-1])
                        nx.add_path(G_direct, [seg_last_false[-1], seg_false_midpoint, seg_last_false[0]])


    G_direct_copy = nx.DiGraph()
    node_start = []
    node_end = []
    for node in G_direct.nodes:
        in_degree = (0 if len({G_direct.in_degree(nbunch=node)}) == 0 else G_direct.in_degree(nbunch=node))
        out_degree = (0 if len({G_direct.out_degree(nbunch=node)}) == 0 else G_direct.out_degree(nbunch=node))
        if in_degree == 0 and out_degree != 0:
            node_start.append(node)
        if in_degree != 0 and out_degree == 0:
            node_end.append(node)
    node_B = []
    for node in g4_copy.nodes:
        if g4_copy.degree[node] > 2 and node not in Disk_node:
            node_B.append(node)
    node_collection = node_B + node_start + node_end
    j0 = 0
    for node in node_collection:
        j0 += 1
        for edge in list(g4_copy.edges(nbunch=node)):
            seg, close_disc = find_seg_new(node, edge[1], close_disc, g4_copy, Disk_node, B_g4_copy)
            if len(seg) == 2:
                continue

            seg_mid1 = seg[int(len(seg)/2)]
            seg_mid2 = seg[-int(len(seg)/2) - 1]
            if seg_mid1 in G_direct.nodes:
                seg_mid = seg_mid1
            if seg_mid2 in G_direct.nodes:
                seg_mid = seg_mid2
            if list(G_direct.in_edges(nbunch=seg_mid))[0][0] == seg[0]:
                path_new = find_path_new(seg)
            if list(G_direct.in_edges(nbunch=seg_mid))[0][0] == seg[-1]:
                seg_1 = seg[::-1]
                path_new = find_path_new(seg_1)
            if path_new[1] not in G_direct_copy:
                nx.add_path(G_direct_copy, path_new)
    return G_direct_copy, G_direct, g4_copy, g4, close_disc