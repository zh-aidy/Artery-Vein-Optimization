import networkx as nx
import math

def find_10x10(i, j, image):
    score = 0
    node_vessel = []
    for m in range(-5, 5):
        for n in range(-5, 5):
            if (image[i + m, j + n, :] != [0, 0, 0]).any():
                score += 1
                node_vessel.append((i + m, j + n))
    return score, node_vessel

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

def find_16x16(i, j, image):
    count = 0
    B_intensive_G_uc_sub = []
    for m in range(-5, 5):
        for n in range(-5, 5):
            if (image[i + m, j + n, :] != [0, 0, 0]).any():
                B_intensive_G_uc_sub.append((i + m, j + n))
                count += 1
    return count, B_intensive_G_uc_sub

def find_seg_passlabel(branchnode, nbr, graph, close_disc):
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

def find_the_path(node, nrb, G_direct_new):
    path = []
    path.append(node)
    path.append(nrb)
    while 1:
        if G_direct_new.degree[nrb] != 2:
            break
        if len(list(G_direct_new.out_edges(nbunch=nrb))) == 0:
            break
        next_point = list(G_direct_new.out_edges(nbunch=nrb))[0][1]
        path.append(next_point)
        if G_direct_new.degree[next_point] != 2:
            break
        else:
            nrb = next_point
    return path

def find_the_path_new(node, nrb, G_direct_new, g4):
    path = []
    path.append(node)
    path.append(nrb)
    while 1:
        if G_direct_new.degree[nrb] != 2 or g4.degree[nrb] != 2:
            break
        if len(list(G_direct_new.out_edges(nbunch=nrb))) == 0:
            break
        next_point = list(G_direct_new.out_edges(nbunch=nrb))[0][1]
        path.append(next_point)
        if G_direct_new.degree[next_point] != 2 or g4.degree[nrb] != 2:
            break
        else:
            nrb = next_point
    return path

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
    if graph.degree[segB[0]] == 1:
        i += 1
    angle = (0.8 ** i) * angle
    return [angle, segA[1]]

def segs_angle_othertwosege(segA, segB):
    angle1 = cal_angle(segA[-1], segA[0], segB[-1])
    angle2 = cal_angle(segA[-1], segA[1], segB[-1])
    angle3 = cal_angle(segA[-1], segB[1], segB[-1])
    angle = 0.8 * angle1 + 0.1 * angle2 + 0.1 * angle3
    return [angle, segA[1]]
def segs_angle_othertwosege1(segA, segB):
    angle1 = cal_angle(segA[0], segA[-1], segB[0])
    angle2 = cal_angle(segA[0], segA[-2], segB[0])
    angle3 = cal_angle(segA[0], segB[-2], segB[0])
    angle = 0.8 * angle1 + 0.1 * angle2 + 0.1 * angle3
    return [angle, segA[1]]

def segs_angle_othertwosege_1(segA, segB):
    angle = cal_angle(segA[1], segA[0], segB[1])
    return [angle, segA[1]]
def segs_angle_othertwosege_2( segA, segB):
    angle1 = cal_angle(segA[2], segA[0], segB[2])
    angle2 = cal_angle(segA[2], segA[1], segB[2])
    angle3 = cal_angle(segA[2], segB[1], segB[2])
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

def find_branchnode_class(G_direct_new, node, g4_not_cut, close_disc, H, Disc_size):
    if G_direct_new.degree[node] == 2:
        return (2)

    if G_direct_new.degree[node] == 3:
        in_degree = (0 if len({G_direct_new.in_degree(nbunch=node)}) == 0 else G_direct_new.in_degree(nbunch=node))
        if in_degree == 1:
            seg = find_seg_passlabel_new(node, list(G_direct_new.in_edges(nbunch=node))[0][0], g4_not_cut, close_disc)
            if seg[-1] in G_direct_new.nodes:
                in_degree_inseg = (0 if len({G_direct_new.in_degree(nbunch=seg[-1])}) == 0 else G_direct_new.in_degree(nbunch=seg[-1]))
                out_degree_inseg = (
                    0 if len({G_direct_new.out_degree(nbunch=seg[-1])}) == 0 else G_direct_new.out_degree(nbunch=seg[-1]))
                if in_degree_inseg == 1 and out_degree_inseg <= 2:
                    if len(seg) < 12 * 2.5:
                        return (3, 1, 0)
                    else:
                        return (3, 1, 1)
                elif in_degree_inseg == 1 and out_degree_inseg > 2:
                    return (3, 1, 3)
                elif in_degree_inseg == 2 and out_degree_inseg == 1:
                    abs_dist = math.sqrt((seg[-1][0] - seg[0][0]) ** 2 + (seg[-1][1] - seg[0][1]) ** 2)
                    if len(seg) < 12 * 2.5 and abs_dist < 20 * 2.5:
                        return (3, 1, 2)
                    else:
                        return (3, 1, 0)



                elif in_degree_inseg == 2 and out_degree_inseg > 1:
                    return (3, 1, 2, 2)
                elif in_degree_inseg == 0:
                    dist = math.sqrt((seg[-1][0] - H[0]) ** 2 + (seg[-1][1] - H[1]) ** 2)
                    if dist < Disc_size * 3:
                        return (3, 3, 1)
                    else:
                        return (3, 0)
                else:
                    return (3)
            else:
                return (3, 1)
        if in_degree == 2:
            return (3, 2)
        else:
            return (3, 3, 2)

    if G_direct_new.degree[node] == 4:
        in_degree = (0 if len({G_direct_new.in_degree(nbunch=node)}) == 0 else G_direct_new.in_degree(nbunch=node))
        if in_degree == 1:
            return (4, 1)
        if in_degree == 2:
            return (4, 2)
        else:
            return (3)
    if G_direct_new.degree[node] == 5:
        return (5)
    else:
        return (3)

def find_twoClasses_cp(G_direct_new, branchnode, g4_not_cut, close_disc):
    in_nrb = []
    out_nrb = []
    [angle_, nrb_] = [0, 0]
    for edge in list(G_direct_new.in_edges(nbunch=branchnode)):
        in_nrb.append(edge[0])
    for edge in list(G_direct_new.out_edges(nbunch=branchnode)):
        out_nrb.append(edge[1])
    nrb = in_nrb + out_nrb
    [nrb_1, nrb_2] = [0, 0]
    for nrb0 in in_nrb:
        for nrb1 in out_nrb:
            seg1 = find_seg_passlabel_new(branchnode, nrb0, g4_not_cut, close_disc)
            seg2 = find_seg_passlabel_new(branchnode, nrb1, g4_not_cut, close_disc)
            [angle2, nrb2] = segs_angle_othertwosege(seg1, seg2)
            if angle2 > angle_:
                [angle_, nrb_] = [angle2, nrb2]
                [nrb_1, nrb_2] = [nrb0, nrb1]
    for node in in_nrb:
        if node != nrb_1:
            nrb_in_2 = node
    for node in out_nrb:
        if node != nrb_2:
            nrb_out_2 = node
    return [(nrb_1, nrb_2), (nrb_in_2, nrb_out_2)]

def find_twoClasses_cp_new(G_direct_new, branchnode, g4_not_cut, close_disc):
    in_nrb = []
    out_nrb = []
    [angle_, nrb_] = [180, 0]
    for edge in list(G_direct_new.in_edges(nbunch=branchnode)):
        in_nrb.append(edge[0])
    for edge in list(G_direct_new.out_edges(nbunch=branchnode)):
        out_nrb.append(edge[1])
    nrb = in_nrb + out_nrb
    [nrb_1, nrb_2] = [0, 0]
    for nrb0 in in_nrb:
        for nrb1 in out_nrb:
            seg1 = find_seg_passlabel_new(branchnode, nrb0, g4_not_cut, close_disc)
            seg2 = find_seg_passlabel_new(branchnode, nrb1, g4_not_cut, close_disc)
            [angle2, nrb2] = segs_angle_othertwosege(seg1, seg2)
            if angle2 < angle_:
                [angle_, nrb_] = [angle2, nrb2]
                [nrb_1, nrb_2] = [nrb0, nrb1]
    for node in in_nrb:
        if node != nrb_1:
            nrb_in_2 = node
    for node in out_nrb:
        if node != nrb_2:
            nrb_out_2 = node
    return [(nrb_1, nrb_out_2), (nrb_in_2, nrb_2), (nrb_1, nrb_2)]

def find_twoClasses_cp_new1(G_direct, node1, g4_not_cut, close_disc):
    nrb_out = list(G_direct.out_edges(nbunch=node1))[0][1]
    seg_out = find_seg_passlabel_new(branchnode=node1, nbr=nrb_out, graph=g4_not_cut, close_disc=close_disc)
    nrb_out_segout_lastnode = []
    for edge in list(G_direct.out_edges(nbunch=seg_out[-1])):
        nrb_out_segout_lastnode.append(edge[1])
    nrb_in = []
    for edge in list(G_direct.in_edges(nbunch=node1)):
        nrb_in.append(edge[0])
    [angle_, nrb_] = [180, 0]
    for nrb1 in nrb_in:
        for nrb2 in nrb_out_segout_lastnode:
            seg_in = find_seg_passlabel_new(branchnode=node1, nbr=nrb1, graph=g4_not_cut, close_disc=close_disc)
            seg_out1 = find_seg_passlabel_new(branchnode=seg_out[-1], nbr=nrb2, graph=g4_not_cut, close_disc=close_disc)
            [angle3, nrb3] = segs_angle_othertwosege(seg_in, seg_out1)
            if angle3 < angle_:
                [angle_, nrb_] = [angle3, nrb3]
                [nrb_1, nrb_2] = [nrb1, nrb2]
    for node in nrb_in:
        if node != nrb_1:
            nrb_in_2 = node
    for node in nrb_out_segout_lastnode:
        if node != nrb_2:
            nrb_out_2 = node
    return [(nrb_1, nrb_out_2), (nrb_in_2, nrb_2), (nrb_1, nrb_2)]

def find_twoClasses_between2node_cp(G_direct_new, node1, g4_not_cut, close_disc):
    nrb_out = []
    for edge in list(G_direct_new.out_edges(nbunch=node1)):
        nrb_out.append(edge[1])
    for edge in list(G_direct_new.in_edges(nbunch=node1)):
        seg_in = find_seg_passlabel_new(node1, edge[0], g4_not_cut, close_disc)
        segIn_midpoint = seg_in[int(len(seg_in)/2)]
    nrb_lastIn = []
    for edge1 in list(G_direct_new.in_edges(nbunch=seg_in[-1])):
        nrb_lastIn.append(edge1[0])
    angle0 = 180
    [nrb_out_final, nrb_lastIn_final] = [0, 0]
    for node in nrb_lastIn:
        for node_ in nrb_out:
            seg_out = find_seg_passlabel_new(node1, node_, g4_not_cut, close_disc)
            seg_last_in = find_seg_passlabel_new(seg_in[-1], node, g4_not_cut, close_disc)
            angle = cal_angle(seg_out[-1], segIn_midpoint, seg_last_in[-1])
            if angle < angle0:
                angle0 = angle
                [nrb_out_final, nrb_lastIn_final] = [node_, node]
    return [nrb_out_final, nrb_lastIn_final]

def find_score(bn_class, G_direct_new, G_direct_v, G_direct_a, node1, g4_not_cut, seg1, nrb_out, close_disc, image_label, H, Disc_size):
    score_v_lastINv = 0
    score_a_lastINa = 0
    if bn_class == (2):
        if len(list(G_direct_v.in_edges(nbunch=node1))) != 0:
            score_v_lastINv = 1
        if len(list(G_direct_a.in_edges(nbunch=node1))) != 0:
            score_a_lastINa = 1
        alpha1 = 0.6
        alpha2 = 0.4
        return alpha1, alpha2, score_v_lastINv, score_a_lastINa

    if bn_class == (3, 1):
        if len(list(G_direct_v.in_edges(nbunch=node1))) != 0:
            score_v_lastINv = 1
        if len(list(G_direct_a.in_edges(nbunch=node1))) != 0:
            score_a_lastINa = 1
        alpha1 = 0.6
        alpha2 = 0.4
        return alpha1, alpha2, score_v_lastINv, score_a_lastINa

    if bn_class == (3, 1, 2, 2):
        nrb_in = list(G_direct_new.in_edges(nbunch=node1))[0][0]
        seg_last = find_seg_passlabel_new(branchnode=node1, nbr=nrb_in, graph=g4_not_cut, close_disc=close_disc)
        [(nrbin_1, nrbout_1), (nrbin_2, nrbout_2), (nrb_1, nrb_2)] = find_twoClasses_cp_new(G_direct_new, seg_last[-1], g4_not_cut, close_disc)
        seg = find_seg_passlabel_new(branchnode=node1, nbr=nrb_out, graph=g4_not_cut, close_disc=close_disc)
        if len(seg_last) < 12 * 2.5:
            if len(list(G_direct_v.in_edges(nbunch=node1))) != 0:
                score_v_lastINv = 1
            if len(list(G_direct_a.in_edges(nbunch=node1))) != 0:
                score_a_lastINa = 1
            if G_direct_new.degree[seg[-1]] == 1:
                alpha1, alpha2 = 0.6, 0.4
            else:
                alpha1 = 0.4
                alpha2 = 0.6
            return alpha1, alpha2, score_v_lastINv, score_a_lastINa

        in_degree_v = (0 if len({G_direct_v.in_degree(nbunch=seg_last[-1])}) == 0 else G_direct_v.in_degree(nbunch=seg_last[-1]))
        in_degree_a = (0 if len({G_direct_a.in_degree(nbunch=seg_last[-1])}) == 0 else G_direct_a.in_degree(
            nbunch=seg_last[-1]))
        out_degree_v = (
            0 if len({G_direct_v.out_degree(nbunch=seg_last[-1])}) == 0 else G_direct_v.out_degree(nbunch=seg_last[-1]))
        out_degree_a = (0 if len({G_direct_a.out_degree(nbunch=seg_last[-1])}) == 0 else G_direct_a.out_degree(
            nbunch=seg_last[-1]))
        if in_degree_v != 0 and in_degree_a != 0:
            if out_degree_v == 0 or out_degree_a == 0:
                if len(list(G_direct_v.in_edges(nbunch=node1))) != 0:
                    score_v_lastINv = 1
                if len(list(G_direct_a.in_edges(nbunch=node1))) != 0:
                    score_a_lastINa = 1
                alpha1 = 0.3
                alpha2 = 0.7
        if nrbout_1 == seg_last[-2]:
            if nrbin_1 in G_direct_v and nrbout_1 in G_direct_v:
                score_v_lastINv = 1
                alpha1 = 0.6
                alpha2 = 0.4
            elif nrbin_1 in G_direct_a and nrbout_1 in G_direct_a:
                score_a_lastINa = 1
                alpha1 = 0.6
                alpha2 = 0.4
            else:
                if len(list(G_direct_v.in_edges(nbunch=node1))) != 0:
                    score_v_lastINv = 1
                if len(list(G_direct_a.in_edges(nbunch=node1))) != 0:
                    score_a_lastINa = 1
                alpha1 = 0.4
                alpha2 = 0.6
        elif nrbout_2 == seg_last[-2]:
            if nrbin_2 in G_direct_v and nrbout_2 in G_direct_v:
                score_v_lastINv = 1
                alpha1 = 0.6
                alpha2 = 0.4
            elif nrbin_2 in G_direct_a and nrbout_2 in G_direct_a:
                score_a_lastINa = 1
                alpha1 = 0.6
                alpha2 = 0.4
            else:
                if len(list(G_direct_v.in_edges(nbunch=node1))) != 0:
                    score_v_lastINv = 1
                if len(list(G_direct_a.in_edges(nbunch=node1))) != 0:
                    score_a_lastINa = 1
                alpha1 = 0.4
                alpha2 = 0.6
        else:
            if len(list(G_direct_v.in_edges(nbunch=node1))) != 0:
                score_v_lastINv = 1
            if len(list(G_direct_a.in_edges(nbunch=node1))) != 0:
                score_a_lastINa = 1
            alpha1 = 0.4
            alpha2 = 0.6
        nrb_inlastlast = []
        for edge in list(G_direct_new.in_edges(nbunch=seg_last[-1])):
            nrb_inlastlast.append(edge[0])
        seg_lastlast1 = find_seg_passlabel_new(branchnode=seg_last[-1], nbr=nrb_inlastlast[0], graph=g4_not_cut, close_disc=close_disc)
        seg_lastlast2 = find_seg_passlabel_new(branchnode=seg_last[-1], nbr=nrb_inlastlast[1], graph=g4_not_cut,
                                               close_disc=close_disc)
        if seg_lastlast1[-1] == seg_lastlast2[-1]:
            alpha1 = 0.3
            alpha2 = 0.4
        return alpha1, alpha2, score_v_lastINv, score_a_lastINa
    if bn_class == (3, 1, 0) or bn_class == (3):
        if len(list(G_direct_v.in_edges(nbunch=node1))) != 0:
            score_v_lastINv = 1
        if len(list(G_direct_a.in_edges(nbunch=node1))) != 0:
            score_a_lastINa = 1
        alpha1 = 0.2
        alpha2 = 0.8
        if bn_class == (3, 1, 0):
            seg_in = find_seg_passlabel_new(node1, (list(G_direct_new.in_edges(nbunch=node1))[0][0]), g4_not_cut, close_disc)
            Seg_in_last = []
            for edge in list(G_direct_new.in_edges(nbunch=seg_in[-1])):
                seg = find_seg_passlabel_new(seg_in[-1], edge[0], g4_not_cut, close_disc)
                Seg_in_last.append(seg)
            if len(Seg_in_last) == 2 and Seg_in_last[0][-1] == Seg_in_last[-1][-1]:
                in_degree_a_last = (0 if len({G_direct_a.in_degree(nbunch=seg_in[-1])}) == 0 else G_direct_a.in_degree(nbunch=seg_in[-1]))
                in_degree_v_last = (0 if len({G_direct_v.in_degree(nbunch=seg_in[-1])}) == 0 else G_direct_v.in_degree(nbunch=seg_in[-1]))
                if (in_degree_a_last == 2 and score_a_lastINa == 1) or (in_degree_v_last == 2 and score_v_lastINv == 1):
                    alpha1, alpha2 = 0.6, 0.4
            seg = find_seg_passlabel_new(node1, nrb_out, g4_not_cut, close_disc)
            if len(Seg_in_last) == 2 and G_direct_new.degree[seg[-1]] == 1:
                alpha1, alpha2 = 0.6, 0.4
        return alpha1, alpha2, score_v_lastINv, score_a_lastINa

    if bn_class == (3, 1, 1):
        seg_in = find_seg_passlabel_new(branchnode=node1, nbr=list(G_direct_new.in_edges(nbunch=node1))[0][0], graph=g4_not_cut,
                                    close_disc=close_disc)
        score_v_lastINv_last = 0
        score_a_lastINa_last = 0
        if len(G_direct_new.in_edges(nbunch=seg_in[-1])) == 1:
            if len(list(G_direct_v.in_edges(nbunch=seg_in[-1]))) != 0:
                score_v_lastINv_last = 1
            if len(list(G_direct_a.in_edges(nbunch=seg_in[-1]))) != 0:
                score_a_lastINa_last = 1
        if len(list(G_direct_v.in_edges(nbunch=node1))) != 0 and score_v_lastINv_last != 0:
            score_v_lastINv = 1 + 0.6 * score_v_lastINv_last
        if len(list(G_direct_a.in_edges(nbunch=node1))) != 0 and score_a_lastINa_last != 0:
            score_a_lastINa = 1 + 0.6 * score_a_lastINa_last
        seg_inlast = find_seg_passlabel_new(branchnode=seg_in[-1], nbr=list(G_direct_new.in_edges(nbunch=seg_in[-1]))[0][0],
                                        graph=g4_not_cut,
                                        close_disc=close_disc)
        nbr_out = []
        for edge in list(G_direct_new.out_edges(nbunch=node1)):
            nbr_out.append(edge[1])
        angle0 = 180
        for node0 in nbr_out:
            seg_out = find_seg_passlabel_new(branchnode=node1, nbr=node0, graph=g4_not_cut, close_disc=close_disc)
            (angle, nrb) = segs_angle_othertwosege(seg_out, seg_in)
            if angle < angle0:
                angle0 = angle
        if G_direct_new.degree[seg_inlast[-1]] == 1 or (score_v_lastINv == 1.6 and seg_inlast[-1] in G_direct_a and
                                                        angle0 < 100) or (score_a_lastINa == 1.6 and seg_inlast[-1] in
                                                                          G_direct_v and angle0 < 100):
            if score_v_lastINv == 1.6:
                score_v_lastINv = 1
            if score_a_lastINa == 1.6:
                score_a_lastINa = 1
        seg_inlast = find_seg_passlabel_new(branchnode=node1,
                                            nbr=list(G_direct_new.in_edges(nbunch=node1))[0][0],
                                            graph=g4_not_cut,
                                            close_disc=close_disc)
        seg_now = find_seg_passlabel_new(branchnode=node1, nbr=nrb_out, graph=g4_not_cut, close_disc=close_disc)
        (angle, nrb) = segs_angle_othertwosege(seg_now, seg_inlast)
        if angle < 90 and len(seg_inlast) < 80 * 2.5:
            if score_v_lastINv == 1.6:
                score_v_lastINv = 1
            if score_a_lastINa == 1.6:
                score_a_lastINa = 1
        if score_v_lastINv == 0 and score_a_lastINa == 0:
            if len(list(G_direct_v.in_edges(nbunch=node1))) != 0:
                score_v_lastINv = 1
            if len(list(G_direct_a.in_edges(nbunch=node1))) != 0:
                score_a_lastINa = 1
        alpha1 = 0.4
        alpha2 = 0.5
        seg = find_seg_passlabel_new(branchnode=node1, nbr=nrb_out, graph=g4_not_cut, close_disc=close_disc)
        score_v, score_a = 0, 0
        i = 0
        j1, j2 = 0, 0
        for node_ in seg:
            i += 1
            if (image_label[node_[0], node_[1], :] == [255, 0, 0]).all():
                score_v += 1
                if i == 1:
                    j1 = 1
            if (image_label[node_[0], node_[1], :] == [0, 0, 255]).all():
                score_a += 1
                if i == 1:
                    j2 = 1
        if j1 == 1:
            score_v = score_v - 2
        if j2 == 1:
            score_a = score_a - 2
        if ((score_v)/(score_v + score_a + 1e-6) > 0.8 and score_v_lastINv == 0) or ((score_a)/(score_v + score_a + 1e-6) > 0.8 and score_a_lastINa == 0):
            X, Y = [], []
            for node_d in G_direct_new.nodes:
                if G_direct_new.degree[node_d] == 1 and G_direct_new.in_degree[node_d] == 1 and node_d not in seg:
                    X.append(node_d[1])
                    Y.append(node_d[0])
            [[(x, y), dist]] = Find_Recent(x=X, y=Y, x_=node1[1], y_=node1[0], k=1)
            if dist < 20 * 2.5:
                nrb_dt = list(G_direct_new.in_edges(nbunch=(y, x)))[0][0]
                seg_dt = find_seg_passlabel_new(branchnode=(y, x), nbr=nrb_dt, graph=g4_not_cut, close_disc=close_disc)
                h0 = node1[0] - y
                w0 = node1[1] - x
                seg_dt_new = []
                for node in seg_dt:
                    seg_dt_new.append((node[0] + h0, node[1] + w0))
                (angle, nrb) = segs_angle_othertwosege(seg, seg_dt_new)
                if angle > 150:
                    if score_a_lastINa == 1.6:
                        score_a_lastINa = 1
                    if score_v_lastINv == 1.6:
                        score_v_lastINv = 1
                    alpha1, alpha2 = 0.3, 0.7
            dist1 = math.sqrt((seg_in[-1][0] - H[0]) ** 2 + (seg_in[-1][1] - H[1]) ** 2)
            if score_a_lastINa + score_v_lastINv == 1.6 and dist1 < Disc_size * 2:
                if score_a_lastINa == 1.6:
                    score_a_lastINa = 1
                if score_v_lastINv == 1.6:
                    score_v_lastINv = 1
        if G_direct_new.out_degree[seg_now[-1]] != 0:
            nrb_nextout = []
            for edge in list(G_direct_new.out_edges(nbunch=seg_now[-1])):
                nrb_nextout.append(edge[1])
            seg_out_next = []
            for node1 in nrb_nextout:
                seg0 = find_seg_passlabel_new(branchnode=seg_now[-1], nbr=node1, graph=g4_not_cut, close_disc=close_disc)
                seg_out_next.append(seg0)
            v, a = 0, 0
            for seg_0 in seg_out_next:
                score_v0, score_a0 = 0, 0
                for node_0 in seg_0:
                    if (image_label[node_0[0], node_0[1], :] == [255, 0, 0]).all():
                        score_v0 += 1
                    if (image_label[node_0[0], node_0[1], :] == [0, 0, 255]).all():
                        score_a0 += 1
                if (score_v0) / (score_v0 + score_a0 + 1e-6) > 0.8:
                    v += 1
                if (score_a0) / (score_v0 + score_a0 + 1e-6) > 0.8:
                    a += 1
            if ((score_v) / (score_v + score_a + 1e-6) > 0.8 and v == len(nrb_nextout) and score_v_lastINv == 0) or ((score_a)/(score_v + score_a + 1e-6) > 0.8 and a == len(nrb_nextout) and score_a_lastINa == 0):
                if score_a_lastINa == 1.6:
                    score_a_lastINa = 1
                if score_v_lastINv == 1.6:
                    score_v_lastINv = 1
        return alpha1, alpha2, score_v_lastINv, score_a_lastINa
    if bn_class == (3, 0):
        if len(list(G_direct_v.in_edges(nbunch=node1))) != 0:
            score_v_lastINv = 1
        if len(list(G_direct_a.in_edges(nbunch=node1))) != 0:
            score_a_lastINa = 1
        alpha1 = 0.5
        alpha2 = 0.4
        seg_inlast = find_seg_passlabel_new(branchnode=node1,
                                            nbr=list(G_direct_new.in_edges(nbunch=node1))[0][0],
                                            graph=g4_not_cut,
                                            close_disc=close_disc)
        seg_now = find_seg_passlabel_new(branchnode=node1, nbr=nrb_out, graph=g4_not_cut, close_disc=close_disc)
        (angle, nrb) = segs_angle_othertwosege(seg_now, seg_inlast)
        if angle < 90:
            alpha1, alpha2 = 0.4, 0.6
        return alpha1, alpha2, score_v_lastINv, score_a_lastINa

    if bn_class == (3, 1, 3):
        if len(list(G_direct_v.in_edges(nbunch=node1))) != 0:
            score_v_lastINv = 1
        if len(list(G_direct_a.in_edges(nbunch=node1))) != 0:
            score_a_lastINa = 1
        nrb_in = list(G_direct_new.in_edges(nbunch=node1))[0][0]
        seg_in = find_seg_passlabel_new(branchnode=node1, nbr=list(G_direct_new.in_edges(nbunch=node1))[0][0],
                                        graph=g4_not_cut,
                                        close_disc=close_disc)
        nrb_last_in = list(G_direct_new.in_edges(nbunch=seg_in[-1]))[0][0]
        i = 0
        if nrb_in in G_direct_a and nrb_last_in in G_direct_a:
            i = 1
            alpha1 = 0.6
            alpha2 = 0.5
        elif nrb_in in G_direct_v and nrb_last_in in G_direct_v:
            i = 2
            alpha1 = 0.6
            alpha2 = 0.5
        else:
            alpha1 = 0.3
            alpha2 = 0.7
        seg_now = find_seg_passlabel_new(branchnode=node1, nbr=nrb_out, graph=g4_not_cut, close_disc=close_disc)
        (angle, nrb) = segs_angle_othertwosege(seg_now, seg_in)
        if angle < 60:
            alpha1 = 0.3
            alpha2 = 0.7
        return alpha1, alpha2, score_v_lastINv, score_a_lastINa

    if bn_class == (3, 3, 1):
        if len(list(G_direct_v.in_edges(nbunch=node1))) != 0:
            score_v_lastINv = 1
        if len(list(G_direct_a.in_edges(nbunch=node1))) != 0:
            score_a_lastINa = 1
        alpha1 = 0.3
        alpha2 = 0.7
        nrbout = []
        for edge in list(G_direct_new.out_edges(nbunch=seg1[-1])):
            nrbout.append(edge[1])
        if nrbout != []:
            score_a = []
            score_v = []
            for node in nrbout:
                seg = find_seg_passlabel_new(branchnode=seg1[-1], nbr=node, graph=g4_not_cut, close_disc=close_disc)
                score_a1, score_v1 = 0, 0
                for node_seg in seg:
                    if (image_label[node_seg[0], node_seg[1], :]== [255, 0, 0]).all():
                        score_v1 += 1
                    if (image_label[node_seg[0], node_seg[1], :] == [0, 0, 255]).all():
                        score_a1 += 1
                score_a.append(score_a1)
                score_v.append(score_v1)
            i1, i2 = 0, 0
            for j in range(0, len(nrbout)):
                if score_a[j] > score_v[j]:
                    i1 += 1
                if score_a[j] < score_v[j]:
                    i2 += 1
            if (score_v_lastINv == 1 and i2 == len(nrbout)) or (score_a_lastINa == 1 and i1 == len(nrbout)):
                alpha1, alpha2 = 0.6, 0.4
        return alpha1, alpha2, score_v_lastINv, score_a_lastINa
    if bn_class == (3, 3, 2):
        if len(list(G_direct_v.in_edges(nbunch=node1))) != 0:
            score_v_lastINv = 1
        if len(list(G_direct_a.in_edges(nbunch=node1))) != 0:
            score_a_lastINa = 1
        alpha1 = 0.3
        alpha2 = 0.7
        return alpha1, alpha2, score_v_lastINv, score_a_lastINa
    if bn_class == (3, 1, 2):
        seg_in = find_seg_passlabel_new(branchnode=node1, nbr=list(G_direct_new.in_edges(nbunch=node1))[0][0],graph=g4_not_cut,close_disc=close_disc)
        [nrb_out_final, nrb_lastIn_final] = find_twoClasses_between2node_cp(G_direct_new, node1, g4_not_cut, close_disc)
        nrb_in_last = []
        for edge in list(G_direct_new.in_edges(nbunch=seg_in[-1])):
            nrb_in_last.append(edge[0])
        seg_in_last1 = find_seg_passlabel_new(branchnode=seg_in[-1], nbr=nrb_in_last[0], graph=g4_not_cut, close_disc=close_disc)
        seg_in_last2 = find_seg_passlabel_new(branchnode=seg_in[-1], nbr=nrb_in_last[1], graph=g4_not_cut,
                                              close_disc=close_disc)
        dist1 = math.sqrt((seg_in_last1[0][0] - seg_in_last1[-1][0]) ** 2 + ((seg_in_last1[0][1] - seg_in_last1[-1][1]) ** 2))
        dist2 = math.sqrt((seg_in_last2[0][0] - seg_in_last2[-1][0]) ** 2 + ((seg_in_last2[0][1] - seg_in_last2[-1][1]) ** 2))
        for node0 in nrb_in_last:
            if node0 != nrb_lastIn_final:
                node_in_another = node0
        Seg_in = [seg_in_last1, seg_in_last2]
        for seg in Seg_in:
            if node_in_another in seg:
                seg_in_another = seg
            else:
                seg_nrb_lastIn_final = seg
        seg_out = find_seg_passlabel_new(branchnode=node1, nbr=nrb_out, graph=g4_not_cut,close_disc=close_disc)
        h0 = seg_in[-1][0] - seg_in[0][0]
        w0 = seg_in[-1][1] - seg_in[0][1]
        seg_in_another_new = []
        for node in seg_in_another:
            seg_in_another_new.append((node[0] + h0, node[1] + w0))
        seg_nrb_lastIn_final_new = []
        for node in seg_nrb_lastIn_final:
            seg_nrb_lastIn_final_new.append((node[0] + h0, node[1] + w0))
        if nrb_lastIn_final in G_direct_a.nodes:
            if nrb_out == nrb_out_final:
                score_v_lastINv = 1
                (angle, nrb) = segs_angle_othertwosege(seg_in_another_new, seg_out)
                if len(seg_in_last1) < 7 * 2.5 and len(seg_in_last2) < 7 * 2.5:
                    alpha1, alpha2 = 0.3, 0.7
                elif dist1 < 10 * 2.5 and dist2 < 10 * 2.5:
                    alpha1, alpha2 = 0.3, 0.7
                elif angle < 100:
                    alpha1, alpha2 = 0.3, 0.7
                elif seg_in_last1[-1] == seg_in_last2[-1]:
                    alpha1, alpha2 = 0.3, 0.7
                else:
                    alpha1, alpha2 = 0.6, 0.4
            else:
                score_a_lastINa = 1
                (angle, nrb) = segs_angle_othertwosege(seg_nrb_lastIn_final_new, seg_out)
                if len(seg_in_last1) < 7 * 2.5 and len(seg_in_last2) < 7 * 2.5:
                    alpha1, alpha2 = 0.3, 0.7
                elif dist1 < 10 * 2.5 and dist2 < 10 * 2.5:
                    alpha1, alpha2 = 0.3, 0.7
                elif angle < 100:
                    alpha1, alpha2 = 0.3, 0.7
                elif seg_in_last1[-1] == seg_in_last2[-1]:
                    alpha1, alpha2 = 0.3, 0.7
                else:
                    alpha1, alpha2 = 0.6, 0.4
        elif nrb_lastIn_final in G_direct_v.nodes:
            if nrb_out == nrb_out_final:
                score_a_lastINa = 1
                (angle, nrb) = segs_angle_othertwosege(seg_in_another_new, seg_out)
                if len(seg_in_last1) < 7 * 2.5 and len(seg_in_last2) < 7 * 2.5:
                    alpha1, alpha2 = 0.3, 0.7
                elif dist1 < 10 * 2.5 and dist2 < 10 * 2.5:
                    alpha1, alpha2 = 0.3, 0.7
                elif angle < 100:
                    alpha1, alpha2 = 0.3, 0.7
                elif seg_in_last1[-1] == seg_in_last2[-1]:
                    alpha1, alpha2 = 0.3, 0.7
                else:
                    alpha1, alpha2 = 0.6, 0.4
            else:
                score_v_lastINv = 1
                (angle, nrb) = segs_angle_othertwosege(seg_nrb_lastIn_final_new, seg_out)
                if len(seg_in_last1) < 7 * 2.5 and len(seg_in_last2) < 7 * 2.5:
                    alpha1, alpha2 = 0.3, 0.7
                elif dist1 < 10 * 2.5 and dist2 < 10 * 2.5:
                    alpha1, alpha2 = 0.3, 0.7
                elif angle < 100:
                    alpha1, alpha2 = 0.3, 0.7
                elif seg_in_last1[-1] == seg_in_last2[-1]:
                    alpha1, alpha2 = 0.3, 0.7
                else:
                    alpha1, alpha2 = 0.6, 0.4
        else:
            alpha1, alpha2 = 0.3, 0.7
            if len(list(G_direct_v.in_edges(nbunch=node1))) != 0:
                score_v_lastINv = 1
            if len(list(G_direct_a.in_edges(nbunch=node1))) != 0:
                score_a_lastINa = 1
        (anglein, nrbin) = segs_angle_othertwosege(seg_in_last1, seg_in_last2)
        nrb_out_list = []
        for edge in list(G_direct_new.out_edges(nbunch=node1)):
            nrb_out_list.append(edge[1])
        seg_out1 = find_seg_passlabel_new(branchnode=node1, nbr=nrb_out_list[0], graph=g4_not_cut, close_disc=close_disc)
        seg_out2 = find_seg_passlabel_new(branchnode=node1, nbr=nrb_out_list[1], graph=g4_not_cut, close_disc=close_disc)
        (angleout, nrbout) = segs_angle_othertwosege(seg_out1, seg_out2)
        if anglein > 105 or angleout > 105:
            alpha1, alpha2 = 0.3, 0.7
        if G_direct_new.in_degree[seg_in_last1[-1]] == 0 or G_direct_new.in_degree[seg_in_last2[-1]] == 0:
            alpha1, alpha2 = 0.3, 0.7
        if seg_out1[-1] == seg_out2[-1]:
            alpha1, alpha2 = 0.3, 0.7
        if dist1 < 20 * 2.5 and dist2 < 20 * 2.5:
            score_v1, score_a1 = 0, 0
            score_v2, score_a2 = 0, 0
            for node in seg_out1:
                if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                    score_v1 += 1
                if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                    score_a1 += 1
            for node in seg_out2:
                if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                    score_v2 += 1
                if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                    score_a2 += 1
            if (score_v1/(score_v1+score_a1+1e-6) > 0.8 and score_a2/(score_a2+score_v2+1e-6)) or ((score_a1/(score_v1+score_a1+1e-6) > 0.8 and score_v2/(score_a2+score_v2+1e-6))):
                alpha1 = 0.4
                alpha2 = 0.6

        return alpha1, alpha2, score_v_lastINv, score_a_lastINa
    if bn_class == (3, 2):
        [angle, nrb] = [0, 0]
        angle = 0
        Seg_in = []
        for edge1 in list(G_direct_new.in_edges(nbunch=node1)):
            seg_v = find_seg_passlabel_new(node1, edge1[0], g4_not_cut, close_disc)
            Seg_in.append(seg_v)
            [angle1, nrb1] = segs_angle_othertwosege(seg_v, seg1)
            if angle1 > angle:
                [angle, nrb] = [angle1, nrb1]
        if nrb in G_direct_v.nodes:
            score_v_lastINv = 1
        if nrb in G_direct_a.nodes:
            score_a_lastINa = 1
        alpha1 = 0.4
        alpha2 = 0.6
        in_degree_a = (0 if len({G_direct_a.in_degree(nbunch=node1)}) == 0 else G_direct_a.in_degree(nbunch=node1))
        in_degree_v = (0 if len({G_direct_v.in_degree(nbunch=node1)}) == 0 else G_direct_v.in_degree(nbunch=node1))
        if Seg_in[0][-1] == Seg_in[-1][-1]:
            nrb_last_in = []
            for edge in list(G_direct_new.in_edges(nbunch=Seg_in[0][-1])):
                nrb_last_in.append(edge[0])
            if len(nrb_last_in) == 1:
                in_degree_a_last = (0 if len({G_direct_a.in_degree(nbunch=Seg_in[0][-1])}) == 0 else G_direct_a.in_degree(nbunch=Seg_in[0][-1]))
                in_degree_v_last = (0 if len({G_direct_v.in_degree(nbunch=Seg_in[0][-1])}) == 0 else G_direct_v.in_degree(nbunch=Seg_in[0][-1]))
                if (in_degree_a == 2 and in_degree_a_last == 1) or (in_degree_v == 2 and in_degree_v_last == 1):
                    alpha1, alpha2 = 0.6, 0.4
        seg = find_seg_passlabel_new(node1, list(G_direct_new.out_edges(nbunch=node1))[0][1], g4_not_cut, close_disc)
        if G_direct_new.degree[seg[-1]] == 1:
            for seg_in in Seg_in:
                if nrb in seg_in:
                    seg_in_last = seg_in
                else:
                    seg_in_another = seg_in
            (i1, j1) = seg_in_last[1]
            (i2, j2) = seg_in_another[1]
            score_another, node_vessel2 = find_10x10(i2, j2, image_label)
            score_last, node_vessel1 = find_10x10(i1, j1, image_label)
            m = 0
            for node in node_vessel1:
                if node in node_vessel2:
                    m += 1
            if score_last > score_another and ((score_last - score_another)/(score_another-m + 1e-6)) > 0.5:
                alpha1, alpha2 = 0.6, 0.4
        return alpha1, alpha2, score_v_lastINv, score_a_lastINa
    if bn_class == (3, 0, 0):
        score_v_lastINv = 0
        score_a_lastINa = 0
        alpha1 = 0
        alpha2 = 1
        return alpha1, alpha2, score_v_lastINv, score_a_lastINa
    if bn_class == (4, 1):
        if len(list(G_direct_v.in_edges(nbunch=node1))) != 0:
            score_v_lastINv = 1
        if len(list(G_direct_a.in_edges(nbunch=node1))) != 0:
            score_a_lastINa = 1
        alpha1 = 0.3
        alpha2 = 0.7
        return alpha1, alpha2, score_v_lastINv, score_a_lastINa
    if bn_class == (4, 2):
        if node1 in G_direct_a and node1 in G_direct_v:
            [(nrbin_1, nrbout_1), (nrbin_2, nrbout_2), (nrb_1, nrb_2)] = find_twoClasses_cp_new(G_direct_new, node1, g4_not_cut, close_disc)
            if nrb_out == nrbout_1:
                if nrbin_1 in G_direct_v.nodes:
                    score_v_lastINv = 1
                if nrbin_1 in G_direct_a.nodes:
                    score_a_lastINa = 1
            if nrb_out == nrbout_2:
                if nrbin_2 in G_direct_v.nodes:
                    score_v_lastINv = 1
                if nrbin_2 in G_direct_a.nodes:
                    score_a_lastINa = 1

            seg1 = find_seg_passlabel_new(node1, nrb_1, g4_not_cut, close_disc)
            seg2 = find_seg_passlabel_new(node1, nrb_2, g4_not_cut, close_disc)
            [angle, nrb_] = segs_angle_othertwosege(seg1, seg2)
            if angle > 130:
                alpha1 = 0.4
                alpha2 = 0.5
            elif angle > 70 and angle <= 130:
                alpha1 = 0.6
                alpha2 = 0.4
            else:
                alpha1 = 0.4
                alpha2 = 0.5
            if alpha1 > alpha2:
                if nrb_out == nrbout_1:
                    seg_out = find_seg_passlabel_new(node1, nrb_out, g4_not_cut, close_disc)
                    seg_in = find_seg_passlabel_new(node1, nrbin_1, g4_not_cut, close_disc)
                    score_a, score_v = 0, 0
                    for node_seg in seg_out:
                        if (image_label[node_seg[0], node_seg[1], :] == [255, 0, 0]).all():
                            score_v += 1
                        if (image_label[node_seg[0], node_seg[1], :] == [0, 0, 255]).all():
                            score_a += 1
                    if (score_a > score_v and nrbin_1 in G_direct_v.nodes) or (score_v > score_a and nrbin_1 in G_direct_a.nodes):
                        dist1 = math.sqrt((seg_out[0][0] - seg_out[-1][0]) ** 2 + (seg_out[0][1] - seg_out[-1][1]) ** 2)
                        dist2 = math.sqrt((seg_in[0][0] - seg_in[-1][0]) ** 2 + (seg_in[0][1] - seg_in[-1][1]) ** 2)
                        Len = min(dist1, dist2)
                        if Len > 20 * 2.5:
                            for node in seg_out:
                                dist = math.sqrt((seg_out[0][0] - node[0]) ** 2 + (seg_out[0][1] - node[1]) ** 2)
                                if dist > 10 * 2.5:
                                    seg_next_0 = node
                                    break
                            for node in seg_in:
                                dist = math.sqrt((seg_in[0][0] - node[0]) ** 2 + (seg_in[0][1] - node[1]) ** 2)
                                if dist > 10 * 2.5:
                                    seg_last_0 = node
                                    break
                            count_out, B_intensive_G_uc_sub = find_16x16(seg_next_0[0], seg_next_0[1], image_label)
                            count_in, B_intensive_G_uc_sub = find_16x16(seg_last_0[0], seg_last_0[1], image_label)
                            if (count_out-count_in)/count_in > 1/3:
                                alpha1, alpha2 = 0.3, 0.7
                if nrb_out == nrbout_2:
                    seg_out = find_seg_passlabel_new(node1, nrb_out, g4_not_cut, close_disc)
                    seg_in = find_seg_passlabel_new(node1, nrbin_2, g4_not_cut, close_disc)
                    score_a, score_v = 0, 0
                    for node_seg in seg_out:
                        if (image_label[node_seg[0], node_seg[1], :] == [255, 0, 0]).all():
                            score_v += 1
                        if (image_label[node_seg[0], node_seg[1], :] == [0, 0, 255]).all():
                            score_a += 1
                    if (score_a > score_v and nrbin_2 in G_direct_v.nodes) or (score_v > score_a and nrbin_2 in G_direct_a.nodes):
                        dist1 = math.sqrt((seg_out[0][0] - seg_out[-1][0]) ** 2 + (seg_out[0][1] - seg_out[-1][1]) ** 2)
                        dist2 = math.sqrt((seg_in[0][0] - seg_in[-1][0]) ** 2 + (seg_in[0][1] - seg_in[-1][1]) ** 2)
                        Len = min(dist1, dist2)
                        if Len > 20 * 2.5:
                            for node in seg_out:
                                dist = math.sqrt((seg_out[0][0] - node[0]) ** 2 + (seg_out[0][1] - node[1]) ** 2)
                                if dist > 10 * 2.5:
                                    seg_next_0 = node
                                    break
                            for node in seg_in:
                                dist = math.sqrt((seg_in[0][0] - node[0]) ** 2 + (seg_in[0][1] - node[1]) ** 2)
                                if dist > 10 * 2.5:
                                    seg_last_0 = node
                                    break
                            count_out, B_intensive_G_uc_sub = find_16x16(seg_next_0[0], seg_next_0[1], image_label)
                            count_in, B_intensive_G_uc_sub = find_16x16(seg_last_0[0], seg_last_0[1], image_label)
                            if (count_out - count_in) / count_in > 1 / 3:
                                alpha1, alpha2 = 0.3, 0.7
        else:
            if node1 in G_direct_a:
                score_a_lastINa = 1
                score_v_lastINv = 0
            if node1 in G_direct_v:
                score_v_lastINv = 1
                score_a_lastINa = 0
            alpha1 = 0.3
            alpha2 = 0.7
        nrb_out_ = []
        for edge in list(G_direct_new.out_edges(nbunch=node1)):
            nrb_out_.append(edge[1])
        seg_out1 = find_seg_passlabel_new(branchnode=node1, nbr=nrb_out_[0], graph=g4_not_cut, close_disc=close_disc)
        seg_out2 = find_seg_passlabel_new(branchnode=node1, nbr=nrb_out_[1], graph=g4_not_cut, close_disc=close_disc)
        if seg_out1[-1] == seg_out2[-1]:
            alpha1 = 0.2
            alpha2 = 0.8
        nrb_in_ = []
        for edge in list(G_direct_new.in_edges(nbunch=node1)):
            nrb_in_.append(edge[0])
        seg_in1 = find_seg_passlabel_new(branchnode=node1, nbr=nrb_in_[0], graph=g4_not_cut, close_disc=close_disc)
        seg_in2 = find_seg_passlabel_new(branchnode=node1, nbr=nrb_in_[1], graph=g4_not_cut, close_disc=close_disc)
        in_degree1 = (
            0 if len({G_direct_new.in_degree(nbunch=seg_in1[-1])}) == 0 else G_direct_new.in_degree(nbunch=seg_in1[-1]))
        in_degree2 = (
            0 if len({G_direct_new.in_degree(nbunch=seg_in2[-1])}) == 0 else G_direct_new.in_degree(nbunch=seg_in2[-1]))
        if in_degree1 == 0 or in_degree2 == 0:
            alpha1 = 0.2
            alpha2 = 0.8
        if seg_in1[-1] == seg_in2[-1]:
            alpha1 = 0.2
            alpha2 = 0.8
        distin1 = math.sqrt((seg_in1[0][0] - seg_in1[-1][0]) ** 2 + ((seg_in1[0][1] - seg_in1[-1][1]) ** 2))
        distin2 = math.sqrt((seg_in2[0][0] - seg_in2[-1][0]) ** 2 + ((seg_in2[0][1] - seg_in2[-1][1]) ** 2))
        if distin1 < 20 * 2.5 and distin2 < 20 * 2.5:
            alpha1 = 0.4
            alpha2 = 0.6
        dist1_2 = math.sqrt((seg_out1[0][0] - seg_out1[-1][0]) ** 2 + ((seg_out1[0][1] - seg_out1[-1][1]) ** 2))
        dist2_2 = math.sqrt((seg_out2[0][0] - seg_out2[-1][0]) ** 2 + ((seg_out2[0][1] - seg_out2[-1][1]) ** 2))
        if dist1_2 < 8 * 2.5 and len(seg_out1) < 5 * 2.5:
            if nrb_out in seg_out2:
                alpha1, alpha2 = 0.4, 0.6
        if dist2_2 < 8 * 2.5 and len(seg_out2) < 5 * 2.5:
            if nrb_out in seg_out1:
                alpha1, alpha2 = 0.4, 0.6
        return alpha1, alpha2, score_v_lastINv, score_a_lastINa
    if bn_class == (5):
        if len(list(G_direct_v.in_edges(nbunch=node1))) != 0:
            for edge1 in list(G_direct_v.in_edges(nbunch=node1)):
                seg_v = find_seg_passlabel_new(node1, edge1[0], g4_not_cut, close_disc)
                [angle1, nrb1] = segs_angle_othertwosege(seg1, seg_v)
                score_v_lastINv += (angle1 - 90) / 90
        if len(list(G_direct_a.in_edges(nbunch=node1))) != 0:
            for edge2 in list(G_direct_a.in_edges(nbunch=node1)):
                seg_a = find_seg_passlabel_new(node1, edge2[0], g4_not_cut, close_disc)
                [angle2, nrb2] = segs_angle_othertwosege(seg1, seg_a)
                score_a_lastINa += (angle2 - 90) / 90
        alpha1 = 0.3
        alpha2 = 0.7
        return alpha1, alpha2, score_v_lastINv, score_a_lastINa

def pass_label(close_disc, G_direct_new, g4, g4_not_cut, image_label, h, H, Disc_size):
    G_direct_v = nx.DiGraph()
    G_direct_a = nx.DiGraph()
    G_direct_v_allseg = nx.DiGraph()
    G_direct_a_allseg = nx.DiGraph()

    j = 0
    N_close_disc_node = len(close_disc)
    miss_1outdegree = []
    G_direct = nx.DiGraph()
    G_direct_allseg = nx.DiGraph()
    BN_G_direct_new = []
    for node in G_direct_new.nodes:
        if G_direct_new.degree[node] > 2:
            BN_G_direct_new.append(node)
    for node0 in close_disc:
        j += 1
        i = 0
        N = len(list(G_direct_new.out_edges(nbunch=node0))[0])
        num_node0nrb = 0
        for edge in list(G_direct_new.out_edges(nbunch=node0)):
            num_node0nrb += 1
            proportion_v_last = 0
            proportion_a_last = 0
            i += 1
            score_v = 0
            score_a = 0
            seg = find_seg_passlabel_new(node0, edge[1], g4_not_cut, close_disc)
            for node_seg in seg:
                if (image_label[node_seg[0], node_seg[1], :] == [255, 0, 0]).all():
                    score_v += 1
                if (image_label[node_seg[0], node_seg[1], :] == [0, 0, 255]).all():
                    score_a += 1
            proportion_v = score_v/len(seg)
            proportion_a = score_a/len(seg)
            path = find_the_path_new(node0, edge[1], G_direct_new, g4_not_cut)

            i0, j0 = 0, 0
            nrb_out = []
            for edge in list(G_direct_new.out_edges(nbunch=seg[-1])):
                nrb_out.append(edge[1])
            if G_direct_new.in_degree[seg[-1]] == 1:
                for node in nrb_out:
                    segout = find_seg_passlabel_new(branchnode=seg[-1], nbr=node, graph=g4_not_cut,
                                                    close_disc=close_disc)
                    score_v_next, score_a_next = 0, 0
                    for node in segout:
                        if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                            score_v_next += 1
                        if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                            score_a_next += 1
                    if score_a_next/(score_v_next+score_a_next+1e-6) > 0.8:
                        i0 += 1
                    if score_v_next/(score_v_next+score_a_next+1e-6) > 0.8:
                        j0 += 1
            if proportion_v > proportion_a:
                if (proportion_a > 0.3 or (len(seg) > 40*2.5 and proportion_a > 0.15))and i0 == len(nrb_out) and i0 != 0:
                    nx.add_path(G_direct_a, path)
                    nx.add_path(G_direct_a_allseg, seg)
                    proportion_v_last = 0
                    proportion_a_last = 1
                else:
                    nx.add_path(G_direct_v, path)
                    nx.add_path(G_direct_v_allseg, seg)
                    proportion_v_last = 1
                    proportion_a_last = 0
            if proportion_v < proportion_a:
                if (proportion_v > 0.3 or (len(seg) > 40*2.5 and proportion_v > 0.15)) and j0 == len(nrb_out) and j0 != 0:
                    nx.add_path(G_direct_v, path)
                    nx.add_path(G_direct_v_allseg, seg)
                    proportion_v_last = 1
                    proportion_a_last = 0
                else:
                    nx.add_path(G_direct_a, path)
                    nx.add_path(G_direct_a_allseg, seg)
                    proportion_v_last = 0
                    proportion_a_last = 1

            nx.add_path(G_direct, path)
            nx.add_path(G_direct_allseg, seg)
            if G_direct_new.degree[path[-1]] == 1:
                if num_node0nrb < G_direct_new.degree[node0]:
                    continue
                if num_node0nrb == G_direct_new.degree[node0]:
                    break
            d = path[-1]
            seg_last = seg
            while 1:
                node = d
                N_node_indegree = (
                    0 if len({G_direct_new.in_degree(nbunch=node)}) == 0 else G_direct_new.in_degree(nbunch=node))
                in_degree = (
                            0 if len({G_direct.in_degree(nbunch=node)}) == 0 else G_direct.in_degree(nbunch=node))
                if N_node_indegree > in_degree:
                    break

                if G_direct_new.degree[node] == 1:
                    break

                out_degree = (
                    0 if len({G_direct_new.out_degree(nbunch=node)}) == 0 else G_direct_new.out_degree(nbunch=node))
                if out_degree == 0:
                    break

                [angle, nrb] = [0, 0]
                for edge in list(G_direct_new.out_edges(nbunch=node)):
                    nrb = edge[1]
                    seg = find_seg_passlabel_new(node, nrb, g4_not_cut, close_disc)
                    seg_0 = [seg[0], seg[1], seg[int(len(seg) / 2)], seg[-2], seg[-1]]
                    seg_last_0 = [seg_last[0], seg_last[1], seg_last[int(len(seg_last) / 2)], seg_last[-2], seg_last[-1]]
                    [angle0, nbr0] = segs_angle_new(seg_0, seg_last_0, g4)
                    if angle0 > angle:
                        [angle, nrb_new] = [angle0, nbr0]
                        seg_next = seg
                score_v_next = 0
                score_a_next = 0
                path_next = find_the_path_new(seg_next[0], seg_next[1], G_direct_new, g4_not_cut)
                for node in seg_next:
                    if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                        score_v_next += 1
                    if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                        score_a_next += 1
                bn_class = find_branchnode_class(G_direct_new, seg_next[0], g4_not_cut, close_disc, H, Disc_size)
                beta1 = 0.4
                beta2 = 0.6
                b1, b2, score_v_lastINv, score_a_lastINa = find_score(bn_class, G_direct_new, G_direct_v,
                                                                      G_direct_a, seg_next[0], g4, seg_next,
                                                                      seg_next[1], close_disc, image_label, H, Disc_size)
                if bn_class in [(2), (3, 1, 0), (3, 1, 1), (3, 1, 3), (3, 1, 2), (3, 1, 2, 2), (3, 3, 1), (3, 3, 2), (3, 0), (3), (3, 1), (3, 2), (4, 1), (4, 2), (5)]:
                    beta1 = b1
                    beta2 = b2
                    proportion_v_last = score_v_lastINv
                    proportion_a_last = score_a_lastINa
                proportion_v_next_0 = score_v_next / (score_v_next + score_a_next + 1e-6)
                proportion_a_next_0 = score_a_next / (score_v_next + score_a_next + 1e-6)
                if len(seg_next) < 10 * 2.5 and proportion_v_next_0 + proportion_a_next_0 < 0.8:
                    proportion_v_next_0 = score_v_next / (score_v_next + score_a_next + 1e-6)
                    proportion_a_next_0 = score_a_next / (score_v_next + score_a_next + 1e-6)
                proportion_v_next = beta1 * proportion_v_last + beta2 * proportion_v_next_0
                proportion_a_next = beta1 * proportion_a_last + beta2 * proportion_a_next_0

                if proportion_v_next > proportion_a_next:
                    nx.add_path(G_direct_v, path_next)
                    nx.add_path(G_direct_v_allseg, seg_next)
                    proportion_v_last = 1
                    proportion_a_last = 0
                if proportion_v_next < proportion_a_next:
                    nx.add_path(G_direct_a, path_next)
                    nx.add_path(G_direct_a_allseg, seg_next)
                    proportion_v_last = 0
                    proportion_a_last = 1
                nx.add_path(G_direct, path_next)
                nx.add_path(G_direct_allseg, seg_next)
                d = path_next[-1]
                seg_last = seg_next
            if i == N:
                break
        if j == N_close_disc_node:
            break

    for node in BN_G_direct_new:
        if node in G_direct.nodes:
            Indegree_G_direct_new = (
                                    0 if len({G_direct_new.in_degree(nbunch=node)}) == 0 else G_direct_new.in_degree(nbunch=node))
            Indegree_G_direct = (
                0 if len({G_direct.in_degree(nbunch=node)}) == 0 else G_direct.in_degree(nbunch=node))
            Outdegree_G_direct_new = (
                0 if len({G_direct_new.out_degree(nbunch=node)}) == 0 else G_direct_new.out_degree(nbunch=node))
            Outdegree_G_direct = (
                0 if len({G_direct.out_degree(nbunch=node)}) == 0 else G_direct.out_degree(nbunch=node))
            if Indegree_G_direct_new == Indegree_G_direct:
                if Outdegree_G_direct_new > Outdegree_G_direct and Indegree_G_direct_new == Indegree_G_direct:
                    miss_1outdegree.append(node)

    BN_G_direct_new = [x for x in BN_G_direct_new if x != (422, 127)]
    while 1:
        Miss_outdgree = miss_1outdegree
        if len(Miss_outdgree) == 0:
            break
        miss_1outdegree = []
        for node1 in Miss_outdgree:
            [angle1, nrb1] = [0, 0]
            score_v_lastINv = 0
            [angle2, nrb2] = [0, 0]
            score_a_lastINa = 0
            [angle1_0, nrb1] = [0, 0]
            score_v_lastINv_0 = 0
            [angle2_0, nrb2] = [0, 0]
            score_a_lastINa_0 = 0

            for edge in list(G_direct_new.out_edges(nbunch=node1)):
                if edge[1] not in G_direct_v.nodes and edge[1] not in G_direct_a.nodes:
                    path = find_the_path_new(node1, edge[1], G_direct_new, g4_not_cut)

                    proportion_v_last = 0
                    proportion_a_last = 0
                    seg1 = find_seg_passlabel_new(node1, edge[1], g4_not_cut, close_disc)
                    bn_class = find_branchnode_class(G_direct_new, seg1[0], g4_not_cut, close_disc, H, Disc_size)
                    score_v_lastINv = 0
                    score_a_lastINa = 0
                    alpha1, alpha2, score_v_lastINv, score_a_lastINa = find_score(bn_class, G_direct_new, G_direct_v, G_direct_a, node1, g4, seg1, edge[1], close_disc, image_label, H, Disc_size)

                    score_v_next = 0
                    score_a_next = 0
                    for node in seg1:
                        if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                            score_v_next += 1
                        if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                            score_a_next += 1
                    proportion_v_next_0 = score_v_next / (score_v_next + score_a_next + 1e-8)
                    proportion_a_next_0 = score_a_next / (score_v_next + score_a_next + 1e-8)
                    if len(seg1) < 10 * 2.5 and proportion_v_next_0 + proportion_a_next_0 < 0.8:
                        proportion_v_next_0 = score_v_next / (score_v_next + score_a_next + 1e-6)
                        proportion_a_next_0 = score_a_next / (score_v_next + score_a_next + 1e-6)
                    proportion_v_next = alpha1 * score_v_lastINv + alpha2 * proportion_v_next_0
                    proportion_a_next = alpha1 * score_a_lastINa + alpha2 * proportion_a_next_0
                    if proportion_v_next > proportion_a_next:
                        nx.add_path(G_direct_v, path)
                        nx.add_path(G_direct_v_allseg, seg1)
                        proportion_v_last = 1
                        proportion_a_last = 0
                    if proportion_v_next < proportion_a_next:
                        nx.add_path(G_direct_a, path)
                        nx.add_path(G_direct_a_allseg, seg1)
                        proportion_v_last = 0
                        proportion_a_last = 1
                    nx.add_path(G_direct, path)
                    nx.add_path(G_direct_allseg, seg1)
            d1 = path[-1]
            seg_last = seg1
            while 1:
                node = d1
                N_node_indegree = (
                    0 if len({G_direct_new.in_degree(nbunch=node)}) == 0 else G_direct_new.in_degree(nbunch=node))
                in_degree = (
                    0 if len({G_direct.in_degree(nbunch=node)}) == 0 else G_direct.in_degree(nbunch=node))
                if N_node_indegree > in_degree:
                    break

                if G_direct_new.degree[node] == 1:
                    break

                out_degree = (
                    0 if len({G_direct_new.out_degree(nbunch=node)}) == 0 else G_direct_new.out_degree(nbunch=node))
                if out_degree == 0:
                    break

                [angle, nrb] = [0, 0]
                for edge in list(G_direct_new.out_edges(nbunch=node)):
                    nrb = edge[1]
                    seg = find_seg_passlabel_new(node, nrb, g4_not_cut, close_disc)
                    seg_0 = [seg[0], seg[1], seg[int(len(seg) / 2)], seg[-2], seg[-1]]
                    seg_last_0 = [seg_last[0], seg_last[1], seg_last[int(len(seg_last) / 2)], seg_last[-2], seg_last[-1]]
                    [angle0, nbr0] = segs_angle_new(seg_0, seg_last_0, g4)
                    if angle0 > angle:
                        [angle, nrb_new] = [angle0, nbr0]
                        seg_next = seg
                score_v_next = 0
                score_a_next = 0
                path_next = find_the_path_new(seg_next[0], seg_next[1], G_direct_new, g4_not_cut)
                for node in seg_next:
                    if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                        score_v_next += 1
                    if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                        score_a_next += 1

                alpha1 = 0.4
                alpha2 = 0.6
                bn_class = find_branchnode_class(G_direct_new, seg_next[0], g4_not_cut, close_disc, H, Disc_size)
                a1, a2, score_v_lastINv, score_a_lastINa = find_score(bn_class, G_direct_new, G_direct_v,
                                                                              G_direct_a, seg_next[0], g4, seg_next, seg_next[1], close_disc, image_label, H, Disc_size)

                if bn_class in [(2), (3, 1, 0), (3, 1, 1), (3, 1, 3), (3, 1, 2), (3, 1, 2, 2), (3, 3, 1), (3, 3, 2), (3, 0), (3),
                                (3, 1), (3, 2), (4, 1), (4, 2), (5)]:
                    alpha1 = a1
                    alpha2 = a2
                    proportion_v_last = score_v_lastINv
                    proportion_a_last = score_a_lastINa
                proportion_v_next_0 = score_v_next / (score_v_next + score_a_next + 1e-6)
                proportion_a_next_0 = score_a_next / (score_v_next + score_a_next + 1e-6)
                if len(seg_next) < 10 * 2.5 and proportion_v_next_0 + proportion_a_next_0 < 0.8:
                    proportion_v_next_0 = score_v_next / (score_v_next + score_a_next + 1e-6)
                    proportion_a_next_0 = score_a_next / (score_v_next + score_a_next + 1e-6)
                proportion_v_next = alpha1 * proportion_v_last + alpha2 * proportion_v_next_0
                proportion_a_next = alpha1 * proportion_a_last + alpha2 * proportion_a_next_0
                if proportion_v_next > proportion_a_next:
                    if bn_class == (4, 2):
                        if proportion_v_last == 1:
                            score_a_next = 0
                            for i in range(0, len(seg_next)):
                                node = seg_next[i]
                                if (image_label[node[0], node[1], :] == [0, 0, 255]).all() and i < 10:
                                    score_a_next += 1
                            if score_a_next > 5 * 2.5:
                                Len = min(len(seg_last), len(seg_next))
                                if Len > 15 * 2.5:
                                    Len = 15 * 2.5
                                seg_next_0 = seg_next[5:Len]
                                seg_last_0 = seg_last[::-1][5:Len]
                                node_seg_next_0 = []
                                for node_0 in seg_next_0:
                                    (i0, j0) = node_0
                                    count, B_intensive_G_uc_sub = find_16x16(i0, j0, image_label)
                                    for node_1 in B_intensive_G_uc_sub:
                                        if node_1 not in node_seg_next_0:
                                            node_seg_next_0.append(node_1)
                                node_seg_last_0 = []
                                for node_0 in seg_last_0:
                                    (i0, j0) = node_0
                                    count, B_intensive_G_uc_sub = find_16x16(i0, j0, image_label)
                                    for node_1 in B_intensive_G_uc_sub:
                                        if node_1 not in node_seg_last_0:
                                            node_seg_last_0.append(node_1)
                                dist_last = 0
                                for i in range(0, len(seg_last_0)-1):
                                    dist_last += math.sqrt((seg_last_0[i][0] - seg_last_0[i+1][0]) ** 2 + (
                                            seg_last_0[i][1] - seg_last_0[i+1][1]) ** 2)
                                dist_next = 0
                                for i in range(0, len(seg_next_0)-1):
                                    dist_next += math.sqrt((seg_next_0[i][0] - seg_next_0[i+1][0]) ** 2 + (
                                            seg_next_0[i][1] - seg_next_0[i+1][1]) ** 2)
                                last = len(node_seg_last_0)/(dist_last + 10)
                                next = len(node_seg_next_0)/(dist_next + 10)
                                if last / (next + 1e-8) < 0.7:
                                    nx.add_path(G_direct_a, path_next)
                                    nx.add_path(G_direct_a_allseg, seg_next)
                                    proportion_a_last = 1
                                    proportion_v_last = 0
                                else:
                                    nx.add_path(G_direct_v, path_next)
                                    nx.add_path(G_direct_v_allseg, seg_next)
                                    proportion_v_last = 1
                                    proportion_a_last = 0
                            else:
                                nx.add_path(G_direct_v, path_next)
                                nx.add_path(G_direct_v_allseg, seg_next)
                                proportion_v_last = 1
                                proportion_a_last = 0
                    else:
                        nx.add_path(G_direct_v, path_next)
                        nx.add_path(G_direct_v_allseg, seg_next)
                        proportion_v_last = 1
                        proportion_a_last = 0
                elif proportion_v_next < proportion_a_next:
                    if bn_class == (4, 2):
                        if proportion_a_last == 1:
                            score_v_next = 0
                            for i in range(0, len(seg_next)):
                                node = seg_next[i]
                                if (image_label[node[0], node[1], :] == [255, 0, 0]).all() and i < 10:
                                    score_v_next += 1
                            if score_v_next > 5 * 2.5:
                                Len = min(len(seg_last), len(seg_next))
                                if Len > 10 * 2.5:
                                    Len = 10 * 2.5
                                seg_next_0 = seg_next[:Len]
                                seg_last_0 = seg_last[::-1][:Len]
                                node_seg_next_0 = []
                                for node_0 in seg_next_0:
                                    (i0, j0) = node_0
                                    count, B_intensive_G_uc_sub = find_16x16(i0, j0, image_label)
                                    for node_1 in B_intensive_G_uc_sub:
                                        if node_1 not in node_seg_next_0:
                                            node_seg_next_0.append(node_1)
                                node_seg_last_0 = []
                                for node_0 in seg_last_0:
                                    (i0, j0) = node_0
                                    count, B_intensive_G_uc_sub = find_16x16(i0, j0, image_label)
                                    for node_1 in B_intensive_G_uc_sub:
                                        if node_1 not in node_seg_last_0:
                                            node_seg_last_0.append(node_1)
                                dist_last = 0
                                for i in range(0, len(seg_last_0) - 1):
                                    dist_last += math.sqrt((seg_last_0[i][0] - seg_last_0[i + 1][0]) ** 2 + (
                                            seg_last_0[i][1] - seg_last_0[i + 1][1]) ** 2)
                                dist_next = 0
                                for i in range(0, len(seg_next_0) - 1):
                                    dist_next += math.sqrt((seg_next_0[i][0] - seg_next_0[i + 1][0]) ** 2 + (
                                            seg_next_0[i][1] - seg_next_0[i + 1][1]) ** 2)
                                last = len(node_seg_last_0) / (dist_last + 10)
                                next = len(node_seg_next_0) / (dist_next + 10)
                                if last / (next + 1e-8) < 0.7:
                                    nx.add_path(G_direct_v, path_next)
                                    nx.add_path(G_direct_v_allseg, seg_next)
                                    proportion_v_last = 1
                                    proportion_a_last = 0
                                else:
                                    nx.add_path(G_direct_a, path_next)
                                    nx.add_path(G_direct_a_allseg, seg_next)
                                    proportion_v_last = 0
                                    proportion_a_last = 1
                        else:
                            nx.add_path(G_direct_a, path_next)
                            nx.add_path(G_direct_a_allseg, seg_next)
                            proportion_v_last = 0
                            proportion_a_last = 1
                    else:
                        nx.add_path(G_direct_a, path_next)
                        nx.add_path(G_direct_a_allseg, seg_next)
                        proportion_v_last = 0
                        proportion_a_last = 1
                else:
                    if proportion_v_next_0 > proportion_a_next_0:
                        nx.add_path(G_direct_v, path_next)
                        nx.add_path(G_direct_v_allseg, seg_next)
                        proportion_v_last = 1
                        proportion_a_last = 0
                    else:
                        nx.add_path(G_direct_a, path_next)
                        nx.add_path(G_direct_a_allseg, seg_next)
                        proportion_v_last = 0
                        proportion_a_last = 1

                nx.add_path(G_direct, path_next)
                nx.add_path(G_direct_allseg, seg_next)
                d1 = path_next[-1]
                seg_last = seg_next

        i1 = 0
        for node in BN_G_direct_new:
            i1 += 1
            if node in G_direct.nodes:
                Indegree_G_direct_new = (
                    0 if len({G_direct_new.in_degree(nbunch=node)}) == 0 else G_direct_new.in_degree(nbunch=node))
                Indegree_G_direct = (
                    0 if len({G_direct.in_degree(nbunch=node)}) == 0 else G_direct.in_degree(nbunch=node))
                Outdegree_G_direct_new = (
                    0 if len({G_direct_new.out_degree(nbunch=node)}) == 0 else G_direct_new.out_degree(nbunch=node))
                Outdegree_G_direct = (
                    0 if len({G_direct.out_degree(nbunch=node)}) == 0 else G_direct.out_degree(nbunch=node))
                if Indegree_G_direct_new == Indegree_G_direct:
                    if Outdegree_G_direct_new > Outdegree_G_direct and Indegree_G_direct_new == Indegree_G_direct:
                        miss_1outdegree.append(node)
            if i1 == len(BN_G_direct_new):
                break

    B4_G_direct = []
    for node in G_direct:
        in_degree = (0 if len({G_direct.in_degree(nbunch=node)}) == 0 else G_direct.in_degree(nbunch=node))
        out_degree = (0 if len({G_direct.out_degree(nbunch=node)}) == 0 else G_direct.out_degree(nbunch=node))
        if in_degree == 1 and out_degree == 3:
            B4_G_direct.append(node)
    for node1 in B4_G_direct:
        nrb_out = []
        for edge in list(G_direct.out_edges(nbunch=node1)):
            nrb_out.append(edge[1])
        for edge in list(G_direct.in_edges(nbunch=node1)):
            nrb_in = edge[0]
            if nrb_in in G_direct_v:
                nrb_out_a = []
                for node2 in nrb_out:
                    if node2 in G_direct_a:
                        nrb_out_a.append(node2)
                if len(nrb_out_a) < 2:
                    break
                segout_anotherclass1 = find_seg_passlabel_new(branchnode=node1, nbr=nrb_out_a[0], graph=g4_not_cut,
                                                          close_disc=close_disc)
                segout_anotherclass1_ = [segout_anotherclass1[0], segout_anotherclass1[1], segout_anotherclass1[int(len(segout_anotherclass1) / 2)], segout_anotherclass1[-2], segout_anotherclass1[-1]]
                segout_anotherclass2 = find_seg_passlabel_new(branchnode=node1, nbr=nrb_out_a[-1], graph=g4_not_cut,
                                                          close_disc=close_disc)
                segout_anotherclass2_ = [segout_anotherclass2[0], segout_anotherclass2[1], segout_anotherclass2[int(len(segout_anotherclass2) / 2)], segout_anotherclass2[-2], segout_anotherclass2[-1]]
                [angle, nrb_in] = segs_angle_othertwosege(segout_anotherclass1_, segout_anotherclass2_)
                if abs(angle - 180) > 30:
                    break
                score0 = -1000
                for node3 in nrb_out_a:
                    seg = find_seg_passlabel_new(branchnode=node1, nbr=node3, graph=g4_not_cut, close_disc=close_disc)
                    seg_lastnode_nrbout = []
                    seg_out = []
                    seg_lastnode_nrbin = []
                    seg_in = []
                    for edge in list(G_direct.out_edges(nbunch=seg[-1])):
                        seg_lastnode_nrbout.append(edge[1])
                        seg_1 = find_seg_passlabel_new(branchnode=seg[-1], nbr=edge[1], graph=g4_not_cut, close_disc=close_disc)
                        seg_out.append([seg_1[0], seg_1[1], seg_1[int(len(seg_1) / 2)], seg_1[-2], seg_1[-1]])
                    for edge in list(G_direct.in_edges(nbunch=seg[-1])):
                        if edge[0] != seg[-2]:
                            seg_lastnode_nrbin.append(edge[0])
                            seg_1 = find_seg_passlabel_new(branchnode=seg[-1], nbr=edge[0], graph=g4_not_cut, close_disc=close_disc)
                            seg_in.append([seg_1[0], seg_1[1], seg_1[int(len(seg_1) / 2)], seg_1[-2], seg_1[-1]])
                    score = 0
                    for seg_2 in seg_in:
                        [angle, nrb_in_] = segs_angle(seg_2, seg)
                        score += 90 - angle
                    for seg_2 in seg_out:
                        [angle, nrb_out_] = segs_angle(seg_2, seg)
                        score += angle - 90
                    if score > score0:
                        seg_true = seg
                        score0 = score
                seg_true_nrb = seg_true[1]
                for node in nrb_out_a:
                    if node != seg_true_nrb:
                        seg_false_nrb = node
                seg_false = find_seg_passlabel_new(branchnode=node1, nbr=seg_false_nrb, graph=g4_not_cut, close_disc=close_disc)
                path_false = find_the_path_new(node=node1, nrb=seg_false_nrb, G_direct_new=G_direct, g4=g4_not_cut)
                path_false_allseg = find_the_path_new(node=node1, nrb=seg_false_nrb, G_direct_new=G_direct_allseg, g4=g4_not_cut)
                if seg_false[-1][0] < h / 2:
                    if seg_false[-1][0] > seg_false[0][0]:
                        seg_false = seg_true
                        path_false = find_the_path_new(node=node1, nrb=seg_false[1], G_direct_new=G_direct,
                                                       g4=g4_not_cut)
                        path_false_allseg = find_the_path_new(node=node1, nrb=seg_false[1],
                                                              G_direct_new=G_direct_allseg, g4=g4_not_cut)
                if seg_false[-1][0] > h / 2:
                    if seg_false[-1][0] < seg_false[0][0]:
                        seg_false = seg_true
                        path_false = find_the_path_new(node=node1, nrb=seg_false[1], G_direct_new=G_direct,
                                                       g4=g4_not_cut)
                        path_false_allseg = find_the_path_new(node=node1, nrb=seg_false[1],
                                                              G_direct_new=G_direct_allseg, g4=g4_not_cut)

                for i in range(0, len(path_false)-1):
                    G_direct.remove_edge(path_false[i], path_false[i+1])
                    G_direct_a.remove_edge(path_false[i], path_false[i+1])
                nx.add_path(G_direct, path_false[::-1])
                nx.add_path(G_direct_a, path_false[::-1])
                for j in range(0, len(path_false_allseg)-1):
                    G_direct_allseg.remove_edge(path_false_allseg[j], path_false_allseg[j+1])
                    G_direct_a_allseg.remove_edge(path_false_allseg[j], path_false_allseg[j+1])
                nx.add_path(G_direct_allseg, path_false_allseg[::-1])
                nx.add_path(G_direct_a_allseg, path_false_allseg[::-1])
            if nrb_in in G_direct_a:
                nrb_out_v = []
                for node2 in nrb_out:
                    if node2 in G_direct_v:
                        nrb_out_v.append(node2)
                if len(nrb_out_v) < 2:
                    break
                segout_anotherclass1 = find_seg_passlabel_new(branchnode=node1, nbr=nrb_out_v[0], graph=g4_not_cut,
                                                          close_disc=close_disc)
                segout_anotherclass1_ = [segout_anotherclass1[0], segout_anotherclass1[1],
                                         segout_anotherclass1[int(len(segout_anotherclass1) / 2)],
                                         segout_anotherclass1[-2], segout_anotherclass1[-1]]
                segout_anotherclass2 = find_seg_passlabel_new(branchnode=node1, nbr=nrb_out_v[-1], graph=g4_not_cut,
                                                          close_disc=close_disc)
                segout_anotherclass2_ = [segout_anotherclass2[0], segout_anotherclass2[1],
                                         segout_anotherclass2[int(len(segout_anotherclass2) / 2)],
                                         segout_anotherclass2[-2], segout_anotherclass2[-1]]
                [angle, nrb_in_] = segs_angle_othertwosege(segout_anotherclass1_, segout_anotherclass2_)
                if abs(angle - 180) > 30:
                    break
                score0 = -1000
                for node3 in nrb_out_v:
                    seg = find_seg_passlabel_new(branchnode=node1, nbr=node3, graph=g4_not_cut, close_disc=close_disc)
                    seg_lastnode_nrbout = []
                    seg_out = []
                    seg_lastnode_nrbin = []
                    seg_in = []
                    for edge in list(G_direct.out_edges(nbunch=seg[-1])):
                        seg_lastnode_nrbout.append(edge[1])
                        seg_1 = find_seg_passlabel_new(branchnode=seg[-1], nbr=edge[1], graph=g4_not_cut, close_disc=close_disc)
                        seg_out.append([seg_1[0], seg_1[1], seg_1[int(len(seg_1) / 2)], seg_1[-2], seg_1[-1]])
                    for edge in list(G_direct.in_edges(nbunch=seg[-1])):
                        if edge[0] != seg[-2]:
                            seg_lastnode_nrbin.append(edge[0])
                            seg_1 = find_seg_passlabel_new(branchnode=seg[-1], nbr=edge[0], graph=g4_not_cut, close_disc=close_disc)
                            seg_in.append([seg_1[0], seg_1[1], seg_1[int(len(seg_1) / 2)], seg_1[-2], seg_1[-1]])
                    score = 0
                    for seg_2 in seg_in:
                        [angle, nrb_in_] = segs_angle(seg_2, seg)
                        score += 90-angle
                    for seg_2 in seg_out:
                        [angle, nrb_out_] = segs_angle(seg_2, seg)
                        score += angle-90
                    if score > score0:
                        seg_true = seg
                        score0 = score
                seg_true_nrb = seg_true[1]
                for node in nrb_out_v:
                    if node != seg_true_nrb:
                        seg_false_nrb = node
                seg_false = find_seg_passlabel_new(branchnode=node1, nbr=seg_false_nrb, graph=g4_not_cut, close_disc=close_disc)
                path_false = find_the_path_new(node=node1, nrb=seg_false_nrb, G_direct_new=G_direct, g4=g4_not_cut)
                path_false_allseg = find_the_path_new(node=node1, nrb=seg_false_nrb, G_direct_new=G_direct_allseg,
                                                      g4=g4_not_cut)
                if seg_false[-1][0] < h/2:
                    if seg_false[-1][0] > seg_false[0][0]:
                        seg_false = seg_true
                        path_false = find_the_path_new(node=node1, nrb=seg_false[1], G_direct_new=G_direct,
                                                       g4=g4_not_cut)
                        path_false_allseg = find_the_path_new(node=node1, nrb=seg_false[1],
                                                              G_direct_new=G_direct_allseg, g4=g4_not_cut)
                if seg_false[-1][0] > h/2:
                    if seg_false[-1][0] < seg_false[0][0]:
                        seg_false = seg_true
                        path_false = find_the_path_new(node=node1, nrb=seg_false[1], G_direct_new=G_direct,
                                                       g4=g4_not_cut)
                        path_false_allseg = find_the_path_new(node=node1, nrb=seg_false[1],
                                                              G_direct_new=G_direct_allseg, g4=g4_not_cut)

                for i in range(0, len(path_false)-1):
                    G_direct.remove_edge(path_false[i], path_false[i+1])
                    G_direct_v.remove_edge(path_false[i], path_false[i+1])
                nx.add_path(G_direct, path_false[::-1])
                nx.add_path(G_direct_v, path_false[::-1])
                for j in range(0, len(path_false_allseg)-1):
                    G_direct_allseg.remove_edge(path_false_allseg[j], path_false_allseg[j+1])
                    G_direct_v_allseg.remove_edge(path_false_allseg[j], path_false_allseg[j+1])
                nx.add_path(G_direct_allseg, path_false_allseg[::-1])
                nx.add_path(G_direct_v_allseg, path_false_allseg[::-1])

    BN_4 = []
    for node in G_direct:
        in_degree = (0 if len({G_direct.in_degree(nbunch=node)}) == 0 else G_direct.in_degree(nbunch=node))
        out_degree = (0 if len({G_direct.out_degree(nbunch=node)}) == 0 else G_direct.out_degree(nbunch=node))
        if in_degree == 2 and out_degree == 2:
            BN_4.append(node)
    for node2 in BN_4:
        out_degree_a = (0 if len({G_direct_a.out_degree(nbunch=node2)}) == 0 else G_direct_a.out_degree(nbunch=node2))
        out_degree_v = (0 if len({G_direct_v.out_degree(nbunch=node2)}) == 0 else G_direct_v.out_degree(nbunch=node2))
        in_degree_a = (0 if len({G_direct_a.in_degree(nbunch=node2)}) == 0 else G_direct_a.in_degree(nbunch=node2))
        in_degree_v = (0 if len({G_direct_v.in_degree(nbunch=node2)}) == 0 else G_direct_v.in_degree(nbunch=node2))
        [(nrbin_1, nrbout_1), (nrbin_2, nrbout_2), (nrb_1, nrb_2)] = find_twoClasses_cp_new(G_direct, node2,
                                                                                            g4_not_cut, close_disc)
        segin1 = find_seg_passlabel_new(branchnode=node2, nbr=nrbin_1, graph=g4_not_cut, close_disc=close_disc)
        segin2 = find_seg_passlabel_new(branchnode=node2, nbr=nrbin_2, graph=g4_not_cut, close_disc=close_disc)
        segout1 = find_seg_passlabel_new(branchnode=node2, nbr=nrbout_1, graph=g4_not_cut, close_disc=close_disc)
        segout2 = find_seg_passlabel_new(branchnode=node2, nbr=nrbout_2, graph=g4_not_cut, close_disc=close_disc)
        (angle__1, nrb__1) = segs_angle_othertwosege(segin1, segin2)
        (angle__2, nrb__1) = segs_angle_othertwosege(segout1, segout2)
        distin1 = math.sqrt((segin1[-1][0] - segin1[0][0]) ** 2 + (segin1[-1][1] - segin1[0][1]) ** 2)
        distin2 = math.sqrt((segin2[-1][0] - segin2[0][0]) ** 2 + (segin2[-1][1] - segin2[0][1]) ** 2)
        distout1 = math.sqrt((segout1[-1][0] - segout1[0][0]) ** 2 + (segout1[-1][1] - segout1[0][1]) ** 2)
        distout2 = math.sqrt((segout2[-1][0] - segout2[0][0]) ** 2 + (segout2[-1][1] - segout2[0][1]) ** 2)
        if out_degree_a == 1 and out_degree_v == 1:
            if segin1[-1] == segin2[-1] or segout1[-1] == segout2[-1]:
                continue

            elif (distin1 < 20 * 2.5 and distin2 < 20 * 2.5) or (distout1 < 20 * 2.5 and distout2 < 20 * 2.5):
                continue
            elif segout1[-1] == segout2[-1]:
                continue
            else:
                if angle__1 > 45 or angle__2 > 45:
                    if nrb_1 in G_direct_a and nrb_2 in G_direct_v:
                        continue
                    elif nrb_1 in G_direct_v and nrb_2 in G_direct_a:
                        continue
                    else:
                        seg1_path_in = find_the_path_new(segin1[-1], segin1[-2], G_direct, g4_not_cut)
                        seg1all_path_in = find_the_path_new(segin1[-1], segin1[-2], G_direct_allseg, g4_not_cut)
                        seg2_path_in = find_the_path_new(segin2[-1], segin2[-2], G_direct, g4_not_cut)
                        seg2all_path_in = find_the_path_new(segin2[-1], segin2[-2], G_direct_allseg, g4_not_cut)
                        if nrbout_1 in G_direct_a and nrbin_1 in G_direct_v:
                            for i in range(0, len(seg1_path_in) - 1):
                                G_direct_v.remove_edge(seg1_path_in[i], seg1_path_in[i + 1])
                            for node1_ in seg1_path_in:
                                G_direct_v.remove_node(node1_)
                            nx.add_path(G_direct_a, seg1_path_in)
                            for i in range(0, len(segin1) - 1):
                                G_direct_v_allseg.remove_edge(seg1all_path_in[i], seg1all_path_in[i + 1])
                            for node1_ in seg1all_path_in:
                                G_direct_v_allseg.remove_node(node1_)
                            nx.add_path(G_direct_a_allseg, seg1all_path_in)
                        if nrbout_1 in G_direct_v and nrbin_1 in G_direct_a:
                            for i in range(0, len(seg1_path_in) - 1):
                                G_direct_a.remove_edge(seg1_path_in[i], seg1_path_in[i + 1])
                            for node1_ in seg1_path_in:
                                G_direct_a.remove_node(node1_)
                            nx.add_path(G_direct_v, seg1_path_in)
                            for i in range(0, len(segin1) - 1):
                                G_direct_a_allseg.remove_edge(seg1all_path_in[i], seg1all_path_in[i + 1])
                            for node1_ in seg1all_path_in:
                                G_direct_a_allseg.remove_node(node1_)
                            nx.add_path(G_direct_v_allseg, seg1all_path_in)
                        if nrbout_2 in G_direct_a and nrbin_2 in G_direct_v:
                            for i in range(0, len(seg2_path_in) - 1):
                                G_direct_v.remove_edge(seg2_path_in[i], seg2_path_in[i + 1])
                            for node1_ in seg2_path_in:
                                G_direct_v.remove_node(node1_)
                            nx.add_path(G_direct_a, seg2_path_in)
                            for i in range(0, len(segin2) - 1):
                                G_direct_v_allseg.remove_edge(seg2all_path_in[i], seg2all_path_in[i + 1])
                            for node1_ in seg2all_path_in:
                                G_direct_v_allseg.remove_node(node1_)
                            nx.add_path(G_direct_a_allseg, seg2all_path_in)
                        if nrbout_2 in G_direct_v and nrbin_2 in G_direct_a:
                            for i in range(0, len(seg2_path_in) - 1):
                                G_direct_a.remove_edge(seg2_path_in[i], seg2_path_in[i + 1])
                            for node1_ in seg2_path_in:
                                G_direct_a.remove_node(node1_)
                            nx.add_path(G_direct_v, seg2_path_in)
                            for i in range(0, len(segin2) - 1):
                                G_direct_a_allseg.remove_edge(seg2all_path_in[i], seg2all_path_in[i + 1])
                            for node1_ in seg2all_path_in:
                                G_direct_a_allseg.remove_node(node1_)
                            nx.add_path(G_direct_v_allseg, seg2all_path_in)
    Bn_oneclass_cross = []
    for node in G_direct:
        in_degree_a = (0 if len({G_direct_a.in_degree(nbunch=node)}) == 0 else G_direct_a.in_degree(nbunch=node))
        in_degree_v = (0 if len({G_direct_v.in_degree(nbunch=node)}) == 0 else G_direct_v.in_degree(nbunch=node))
        if in_degree_a == 2 or in_degree_v == 2:
            Bn_oneclass_cross.append(node)
    for node1 in Bn_oneclass_cross:
        if node1 == (230, 219):
            continue
        in_degree_a = (0 if len({G_direct_a.in_degree(nbunch=node1)}) == 0 else G_direct_a.in_degree(nbunch=node1))
        in_degree_v = (0 if len({G_direct_v.in_degree(nbunch=node1)}) == 0 else G_direct_v.in_degree(nbunch=node1))
        bn_class = find_branchnode_class(G_direct_new=G_direct, node=node1, g4_not_cut=g4_not_cut, close_disc=close_disc, H=H, Disc_size=Disc_size)
        if in_degree_a == 2 or in_degree_v == 2:
            if g4.degree[node1] == 3:
                nrb_in__ = []
                for edge in list(G_direct_allseg.in_edges(nbunch=node1)):
                    nrb_in__.append(edge[0])
                seg1_ = find_seg_passlabel_new(branchnode=node1, nbr=nrb_in__[0], graph=g4_not_cut, close_disc=close_disc)
                seg2_ = find_seg_passlabel_new(branchnode=node1, nbr=nrb_in__[1], graph=g4_not_cut,close_disc=close_disc)
                if seg1_[-1] == seg2_[-1]:
                    score_segin1_v = 0
                    score_segin1_a = 0
                    score_segin2_v = 0
                    score_segin2_a = 0
                    for node in seg1_:
                        if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                            score_segin1_v += 1
                        if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                            score_segin1_a += 1
                    for node in seg2_:
                        if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                            score_segin2_v += 1
                        if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                            score_segin2_a += 1
                    p_a1 = score_segin1_a/(score_segin1_v + score_segin1_a + 1e-6)
                    p_v1 = score_segin1_v/(score_segin1_v + score_segin1_a + 1e-6)
                    p_a2 = score_segin2_a/(score_segin2_v + score_segin2_a + 1e-6)
                    p_v2 = score_segin2_v / (score_segin2_v + score_segin2_a + 1e-6)
                    if (p_a1 > 0.8 and p_a2 > 0.8) or (p_v1 > 0.8 and p_v2 > 0.8):
                        continue

            if bn_class == (4, 2):
                [(nrbin_1, nrbout_1), (nrbin_2, nrbout_2), (nrb_1, nrb_2)] = find_twoClasses_cp_new(G_direct, node1, g4,
                                                                                                close_disc)
                segin_1 = find_seg_passlabel_new(branchnode=node1, nbr=nrbin_1, graph=g4_not_cut, close_disc=close_disc)
                segout_1 = find_seg_passlabel_new(branchnode=node1, nbr=nrbout_1, graph=g4_not_cut, close_disc=close_disc)
                segin_2 = find_seg_passlabel_new(branchnode=node1, nbr=nrbin_2, graph=g4_not_cut, close_disc=close_disc)
                segout_2 = find_seg_passlabel_new(branchnode=node1, nbr=nrbout_2, graph=g4_not_cut, close_disc=close_disc)
                [angle, nrb] = [0, 0]
                for edge in list(G_direct.in_edges(nbunch=segin_1[-1])):
                    segin_1_last1 = find_seg_passlabel_new(branchnode=segin_1[-1], nbr=edge[0], graph=g4_not_cut,
                                                       close_disc=close_disc)
                    segin_1_last1_ = [segin_1_last1[0], segin_1_last1[1], segin_1_last1[int(len(segin_1_last1) / 2)], segin_1_last1[-2], segin_1_last1[-1]]
                    segin_1_ = [segin_1[0], segin_1[1], segin_1[int(len(segin_1) / 2)], segin_1[-2], segin_1[-1]]
                    [angle0, nbr0] = segs_angle_new(segin_1_last1_, segin_1_, g4)
                    if angle0 > angle:
                        [angle, nrb_new] = [angle0, nbr0]
                segin_1_lastIn = find_seg_passlabel_new(branchnode=segin_1[-1], nbr=nrb_new, graph=g4_not_cut,
                                                    close_disc=close_disc)
                [angle1, nrb1] = [0, 0]
                for edge in list(G_direct.in_edges(nbunch=segin_2[-1])):
                    segin_2_last1 = find_seg_passlabel_new(branchnode=segin_2[-1], nbr=edge[0], graph=g4_not_cut,
                                                       close_disc=close_disc)
                    segin_2_last1_ = [segin_2_last1[0], segin_2_last1[1], segin_2_last1[int(len(segin_2_last1) / 2)], segin_2_last1[-2], segin_2_last1[-1]]
                    segin_2_ = [segin_2[0], segin_2[1], segin_2[int(len(segin_2) / 2)], segin_2[-2], segin_2[-1]]
                    [angle0, nbr0] = segs_angle_new(segin_2_last1_, segin_2_, g4)
                    if angle0 > angle1:
                        [angle1, nrb_new] = [angle0, nbr0]
                segin_2_lastIn = find_seg_passlabel_new(branchnode=segin_2[-1], nbr=nrb_new, graph=g4_not_cut,
                                                    close_disc=close_disc)
                score_segin1_v = 0
                score_segin1_a = 0
                segin_1_lastIn_v = 0
                segin_1_lastIn_a = 0
                score_segout1_v = 0
                score_segout1_a = 0
                segin_2_lastIn_v = 0
                segin_2_lastIn_a = 0
                score_segin2_v = 0
                score_segin2_a = 0
                score_segout2_v = 0
                score_segout2_a = 0
                for node in segin_1:
                    if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                        score_segin1_v += 1
                    if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                        score_segin1_a += 1
                proportion_v_segin1 = score_segin1_v / len(segin_1)
                proportion_a_segin1 = score_segin1_a / len(segin_1)
                for node in segin_1_lastIn:
                    if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                        segin_1_lastIn_v += 1
                    if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                        segin_1_lastIn_a += 1
                proportion_v_segin_1_lastIn = segin_1_lastIn_v / len(segin_1_lastIn)
                proportion_a_segin_1_lastIn = segin_1_lastIn_a / len(segin_1_lastIn)
                for node in segout_1:
                    if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                        score_segout1_v += 1
                    if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                        score_segout1_a += 1
                proportion_v_segout1 = score_segout1_v / len(segout_1)
                proportion_a_segout1 = score_segout1_a / len(segout_1)
                for node in segin_2:
                    if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                        score_segin2_v += 1
                    if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                        score_segin2_a += 1
                proportion_v_segin2 = score_segin2_v / len(segin_2)
                proportion_a_segin2 = score_segin2_a / len(segin_2)
                for node in segin_2_lastIn:
                    if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                        segin_2_lastIn_v += 1
                    if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                        segin_2_lastIn_a += 1
                proportion_v_segin_2_lastIn = segin_2_lastIn_v / len(segin_2_lastIn)
                proportion_a_segin_2_lastIn = segin_2_lastIn_a / len(segin_2_lastIn)
                for node in segout_2:
                    if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                        score_segout2_v += 1
                    if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                        score_segout2_a += 1
                proportion_v_segout2 = score_segout2_v / len(segout_2)
                proportion_a_segout2 = score_segout2_a / len(segout_2)
                seg1_v = proportion_v_segin1 + proportion_v_segout1 + proportion_v_segin_1_lastIn
                seg1_a = proportion_a_segin1 + proportion_a_segout1 + proportion_a_segin_1_lastIn
                seg2_v = proportion_v_segin2 + proportion_v_segout2 + proportion_v_segin_2_lastIn
                seg2_a = proportion_a_segin2 + proportion_a_segout2 + proportion_a_segin_2_lastIn
                if in_degree_v == 2:
                    if seg1_v > seg2_v:
                        seg2_path_in = find_the_path_new(segin_2[-1], segin_2[-2], G_direct, g4_not_cut)
                        seg2_path_out = find_the_path_new(node1, nrbout_2, G_direct, g4_not_cut)
                        seg2all_path_in = find_the_path_new(segin_2[-1], segin_2[-2], G_direct_allseg, g4_not_cut)
                        seg2all_path_out = find_the_path_new(node1, nrbout_2, G_direct_allseg, g4_not_cut)
                        if seg2_path_in[-2] in G_direct_v:
                            for node1_ in seg2_path_in[1: -1]:
                                G_direct_v.remove_node(node1_)
                            nx.add_path(G_direct_a, seg2_path_in)
                            for node1_ in seg2all_path_in[1: -1]:
                                G_direct_v_allseg.remove_node(node1_)
                            nx.add_path(G_direct_a_allseg, seg2all_path_in)
                        seg1_path_out = find_the_path_new(node1, nrbout_1, G_direct, g4_not_cut)
                        if seg1_path_out[1] in G_direct_a and seg1_path_out[-1] == seg2_path_out[-1]:
                            continue
                        seg1_path_in = find_the_path_new(segin_1[-1], segin_1[-2], G_direct, g4_not_cut)
                        if seg1_path_in[0] == seg2_path_in[0] or seg1_path_out[-1] == seg2_path_out[-1]:
                            continue

                        seg1all_path_in = find_the_path_new(segin_1[-1], segin_1[-2], G_direct_allseg, g4_not_cut)
                        seg1all_path_out = find_the_path_new(node1, nrbout_1, G_direct_allseg, g4_not_cut)
                        seg2all_path_in_ = [seg2all_path_in[0], seg2all_path_in[1], seg2all_path_in[int(len(seg2all_path_in) / 2)], seg2all_path_in[-2],seg2all_path_in[-1]]
                        seg2all_path_out_ = [seg2all_path_out[0], seg2all_path_out[1], seg2all_path_out[int(len(seg2all_path_out) / 2)], seg2all_path_out[-2],seg2all_path_out[-1]]
                        seg1all_path_in_ = [seg1all_path_in[0], seg1all_path_in[1], seg1all_path_in[int(len(seg1all_path_in) / 2)], seg1all_path_in[-2], seg1all_path_in[-1]]
                        seg1all_path_out_ = [seg1all_path_out[0], seg1all_path_out[1],seg1all_path_out[int(len(seg1all_path_out) / 2)], seg1all_path_out[-2],seg1all_path_out[-1]]
                        (anglein, nrbin) = segs_angle_othertwosege1(seg2all_path_in_, seg1all_path_in_)
                        (angleout, nrbout) = segs_angle_othertwosege(seg2all_path_out_, seg1all_path_out_)
                        if angleout > 20 or anglein > 20:
                            if seg2_path_out[2] in G_direct_v:
                                for node1_ in seg2_path_out[1: -1]:
                                    G_direct_v.remove_node(node1_)
                                nx.add_path(G_direct_a, seg2_path_out)
                                for node1_ in seg2all_path_out[1: -1]:
                                    G_direct_v_allseg.remove_node(node1_)
                                nx.add_path(G_direct_a_allseg, seg2all_path_out)
                            if seg1_path_out[1] in G_direct_a:
                                a1, a2, score_v_lastINv, score_a_lastINa = find_score(bn_class, G_direct_new,
                                                                                      G_direct_v,
                                                                                      G_direct_a, seg1_path_out[0], g4,
                                                                                      seg1all_path_out,
                                                                                      seg1all_path_out[1], close_disc, image_label, H, Disc_size)
                                if a1 > a2:
                                    for node1_ in seg1_path_out[1: -1]:
                                        G_direct_a.remove_node(node1_)
                                    nx.add_path(G_direct_v, seg1_path_out)
                                    for node1_ in seg1all_path_out[1: -1]:
                                        G_direct_a_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_v_allseg, seg1all_path_out)
                                if G_direct_a.out_degree[seg1_path_out[-1]] != 0:
                                    nrb_out_a_ = list(G_direct_a.out_edges(nbunch=seg1_path_out[-1]))[0][1]
                                    seg_a_ = find_seg_passlabel_new(branchnode=seg1_path_out[-1], nbr=nrb_out_a_,
                                                                    graph=g4_not_cut, close_disc=close_disc)
                                    score_v_next = 0
                                    score_a_next = 0
                                    for node in seg_a_:
                                        if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                                            score_v_next += 1
                                        if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                                            score_a_next += 1
                                    seg_a_all_path_out = find_the_path_new(seg1_path_out[-1], nrb_out_a_,
                                                                           G_direct_allseg, g4_not_cut)
                                    seg_a_path_out = find_the_path_new(seg1_path_out[-1], nrb_out_a_, G_direct,
                                                                       g4_not_cut)
                                    if score_v_next > score_a_next:
                                        for node1_ in seg_a_path_out[1: -1]:
                                            G_direct_a.remove_node(node1_)
                                        nx.add_path(G_direct_v, seg_a_path_out)
                                        for node1_ in seg_a_all_path_out[1: -1]:
                                            G_direct_a_allseg.remove_node(node1_)
                                        nx.add_path(G_direct_v_allseg, seg_a_all_path_out)
                                if G_direct.out_degree[seg1_path_out[-1]] == G_direct_v.out_degree[
                                    seg1_path_out[-1]]:
                                    if seg1_path_out[1] in G_direct_a:
                                        for node1_ in seg1_path_out[1: -1]:
                                            G_direct_a.remove_node(node1_)
                                        nx.add_path(G_direct_v, seg1_path_out)
                                        for node1_ in seg1all_path_out[1: -1]:
                                            G_direct_a_allseg.remove_node(node1_)
                                        nx.add_path(G_direct_v_allseg, seg1all_path_out)
                        else:
                            if seg1_path_out[1] in G_direct_v:
                                for node1_ in seg1_path_out[1: -1]:
                                    G_direct_v.remove_node(node1_)
                                nx.add_path(G_direct_a, seg1_path_out)
                                for node1_ in seg1all_path_out[1: -1]:
                                    G_direct_v_allseg.remove_node(node1_)
                                nx.add_path(G_direct_a_allseg, seg1all_path_out)
                            if seg2_path_out[2] in G_direct_v:
                                continue
                            else:
                                for node1_ in seg2_path_out[1: -1]:
                                    G_direct_a.remove_node(node1_)
                                nx.add_path(G_direct_v, seg2_path_out)
                                for node1_ in seg2all_path_out[1: -1]:
                                    G_direct_a_allseg.remove_node(node1_)
                                nx.add_path(G_direct_v_allseg, seg2all_path_out)

                    else:
                        seg1_path_in = find_the_path_new(segin_1[-1], segin_1[-2], G_direct, g4_not_cut)
                        seg1_path_out = find_the_path_new(node1, nrbout_1, G_direct, g4_not_cut)
                        seg1all_path_in = find_the_path_new(segin_1[-1], segin_1[-2], G_direct_allseg, g4_not_cut)
                        seg1all_path_out = find_the_path_new(node1, nrbout_1, G_direct_allseg, g4_not_cut)
                        if seg1_path_in[-2] in G_direct_v:
                            for node1_ in seg1_path_in[1: -1]:
                                G_direct_v.remove_node(node1_)
                            nx.add_path(G_direct_a, seg1_path_in)
                            for node1_ in seg1all_path_in[1: -1]:
                                G_direct_v_allseg.remove_node(node1_)
                            nx.add_path(G_direct_a_allseg, seg1all_path_in)
                        seg2_path_out = find_the_path_new(node1, nrbout_2, G_direct, g4_not_cut)
                        if seg2_path_out[1] in G_direct_a and seg2_path_out[-1] == seg1_path_out[-1]:
                            continue
                        seg2_path_in = find_the_path_new(segin_2[-1], segin_2[-2], G_direct, g4_not_cut)
                        if seg1_path_in[0] == seg2_path_in[0] or seg1_path_out[-1] == seg2_path_out[-1]:
                            continue

                        seg2all_path_in = find_the_path_new(segin_2[-1], segin_2[-2], G_direct_allseg, g4_not_cut)
                        seg2all_path_out = find_the_path_new(node1, nrbout_2, G_direct_allseg, g4_not_cut)
                        seg2all_path_in_ = [seg2all_path_in[0], seg2all_path_in[1], seg2all_path_in[int(len(seg2all_path_in) / 2)], seg2all_path_in[-2], seg2all_path_in[-1]]
                        seg2all_path_out_ = [seg2all_path_out[0], seg2all_path_out[1], seg2all_path_out[int(len(seg2all_path_out) / 2)], seg2all_path_out[-2], seg2all_path_out[-1]]
                        seg1all_path_in_ = [seg1all_path_in[0], seg1all_path_in[1], seg1all_path_in[int(len(seg1all_path_in) / 2)], seg1all_path_in[-2], seg1all_path_in[-1]]
                        seg1all_path_out_ = [seg1all_path_out[0], seg1all_path_out[1],
                                             seg1all_path_out[int(len(seg1all_path_out) / 2)], seg1all_path_out[-2],
                                             seg1all_path_out[-1]]
                        (anglein, nrbin) = segs_angle_othertwosege1(seg2all_path_in_, seg1all_path_in_)
                        (angleout, nrbout) = segs_angle_othertwosege(seg2all_path_out_, seg1all_path_out_)
                        if angleout > 20 or anglein > 20:
                            if seg1_path_out[2] in G_direct_v:
                                for node1_ in seg1_path_out[1: -1]:
                                    G_direct_v.remove_node(node1_)
                                nx.add_path(G_direct_a, seg1_path_out)
                                for node1_ in seg1all_path_out[1: -1]:
                                    G_direct_v_allseg.remove_node(node1_)
                                nx.add_path(G_direct_a_allseg, seg1all_path_out)
                            seg2all_path_out = find_the_path_new(node1, nrbout_2, G_direct_allseg, g4_not_cut)
                            if seg2_path_out[1] in G_direct_a:
                                a1, a2, score_v_lastINv, score_a_lastINa = find_score(bn_class, G_direct_new, G_direct_v,
                                                                                      G_direct_a, seg1_path_out[0], g4, seg2all_path_out,
                                                                                      seg2all_path_out[1], close_disc, image_label, H,Disc_size)
                                if a1 > a2:
                                    for node1_ in seg2_path_out[1: -1]:
                                        G_direct_a.remove_node(node1_)
                                    nx.add_path(G_direct_v, seg2_path_out)
                                    for node1_ in seg2all_path_out[1: -1]:
                                        G_direct_a_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_v_allseg, seg2all_path_out)
                                if G_direct_a.out_degree[seg2_path_out[-1]] != 0:
                                    nrb_out_a_ = list(G_direct_a.out_edges(nbunch=seg2_path_out[-1]))[0][1]
                                    seg_a_ = find_seg_passlabel_new(branchnode=seg2_path_out[-1], nbr=nrb_out_a_, graph=g4_not_cut, close_disc=close_disc)
                                    score_v_next = 0
                                    score_a_next = 0
                                    for node in seg_a_:
                                        if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                                            score_v_next += 1
                                        if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                                            score_a_next += 1
                                    seg_a_all_path_out = find_the_path_new(seg2_path_out[-1], nrb_out_a_, G_direct_allseg, g4_not_cut)
                                    seg_a_path_out = find_the_path_new(seg2_path_out[-1], nrb_out_a_, G_direct, g4_not_cut)
                                    if score_v_next > score_a_next:
                                        for node1_ in seg_a_path_out[1: -1]:
                                            G_direct_a.remove_node(node1_)
                                        nx.add_path(G_direct_v, seg_a_path_out)
                                        for node1_ in seg_a_all_path_out[1: -1]:
                                            G_direct_a_allseg.remove_node(node1_)
                                        nx.add_path(G_direct_v_allseg, seg_a_all_path_out)
                                if G_direct.out_degree[seg2_path_out[-1]] == G_direct_v.out_degree[seg2_path_out[-1]]:
                                    if seg2_path_out[1] in G_direct_a:
                                        for node1_ in seg2_path_out[1: -1]:
                                            G_direct_a.remove_node(node1_)
                                        nx.add_path(G_direct_v, seg2_path_out)
                                        for node1_ in seg2all_path_out[1: -1]:
                                            G_direct_a_allseg.remove_node(node1_)
                                        nx.add_path(G_direct_v_allseg, seg2all_path_out)
                        else:
                            if seg2_path_out[1] in G_direct_v:
                                for node1_ in seg2_path_out[1: -1]:
                                    G_direct_v.remove_node(node1_)
                                nx.add_path(G_direct_a, seg2_path_out)
                                for node1_ in seg2all_path_out[1: -1]:
                                    G_direct_v_allseg.remove_node(node1_)
                                nx.add_path(G_direct_a_allseg, seg2all_path_out)
                            if seg1_path_out[2] in G_direct_v:
                                continue
                            else:
                                for node1_ in seg1_path_out[1: -1]:
                                    G_direct_a.remove_node(node1_)
                                nx.add_path(G_direct_v, seg1_path_out)
                                for node1_ in seg1all_path_out[1: -1]:
                                    G_direct_a_allseg.remove_node(node1_)
                                nx.add_path(G_direct_v_allseg, seg1all_path_out)

                if in_degree_a == 2:
                    if seg1_a > seg2_a:
                        seg2_path_in = find_the_path_new(segin_2[-1], segin_2[-2], G_direct, g4_not_cut)
                        seg2_path_out = find_the_path_new(node1, nrbout_2, G_direct, g4_not_cut)
                        seg2all_path_in = find_the_path_new(segin_2[-1], segin_2[-2], G_direct_allseg, g4_not_cut)
                        seg2all_path_out = find_the_path_new(node1, nrbout_2, G_direct_allseg, g4_not_cut)
                        if seg2_path_in[-2] in G_direct_a:
                            for node1_ in seg2_path_in[1: -1]:
                                G_direct_a.remove_node(node1_)
                            nx.add_path(G_direct_v, seg2_path_in)
                            for node1_ in seg2all_path_in[1: -1]:
                                G_direct_a_allseg.remove_node(node1_)
                            nx.add_path(G_direct_v_allseg, seg2all_path_in)
                        seg1_path_out = find_the_path_new(node1, nrbout_1, G_direct, g4_not_cut)
                        if seg1_path_out[1] in G_direct_a and seg1_path_out[-1] == seg2_path_out[-1]:
                            continue
                        seg1_path_in = find_the_path_new(segin_1[-1], segin_1[-2], G_direct, g4_not_cut)
                        if seg1_path_in[0] == seg2_path_in[0] or seg1_path_out[-1] == seg2_path_out[-1]:
                            continue

                        seg1all_path_in = find_the_path_new(segin_1[-1], segin_1[-2], G_direct_allseg, g4_not_cut)
                        seg1all_path_out = find_the_path_new(node1, nrbout_1, G_direct_allseg, g4_not_cut)
                        seg2all_path_in_ = [seg2all_path_in[0], seg2all_path_in[1],
                                            seg2all_path_in[int(len(seg2all_path_in) / 2)], seg2all_path_in[-2],
                                            seg2all_path_in[-1]]
                        seg2all_path_out_ = [seg2all_path_out[0], seg2all_path_out[1],
                                             seg2all_path_out[int(len(seg2all_path_out) / 2)], seg2all_path_out[-2],
                                             seg2all_path_out[-1]]
                        seg1all_path_in_ = [seg1all_path_in[0], seg1all_path_in[1],
                                            seg1all_path_in[int(len(seg1all_path_in) / 2)], seg1all_path_in[-2],
                                            seg1all_path_in[-1]]
                        seg1all_path_out_ = [seg1all_path_out[0], seg1all_path_out[1],
                                             seg1all_path_out[int(len(seg1all_path_out) / 2)], seg1all_path_out[-2],
                                             seg1all_path_out[-1]]
                        (anglein, nrbin) = segs_angle_othertwosege1(seg2all_path_in_, seg1all_path_in_)
                        (angleout, nrbout) = segs_angle_othertwosege(seg2all_path_out_, seg1all_path_out_)
                        if angleout > 20 or anglein > 20:
                            if seg2_path_out[2] in G_direct_a:
                                for node1_ in seg2_path_out[1: -1]:
                                    G_direct_a.remove_node(node1_)
                                nx.add_path(G_direct_v, seg2_path_out)
                                for node1_ in seg2all_path_out[1: -1]:
                                    G_direct_a_allseg.remove_node(node1_)
                                nx.add_path(G_direct_v_allseg, seg2all_path_out)
                            if seg1_path_out[1] in G_direct_v:
                                a1, a2, score_v_lastINv, score_a_lastINa = find_score(bn_class, G_direct_new,
                                                                                      G_direct_v,
                                                                                      G_direct_a, seg1_path_out[0], g4,
                                                                                      seg1all_path_out,
                                                                                      seg1all_path_out[1], close_disc, image_label, H, Disc_size)
                                if a1 > a2:
                                    for node1_ in seg1_path_out[1: -1]:
                                        G_direct_v.remove_node(node1_)
                                    nx.add_path(G_direct_a, seg1_path_out)
                                    for node1_ in seg1all_path_out[1: -1]:
                                        G_direct_v_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_a_allseg, seg1all_path_out)
                                if G_direct_v.out_degree[seg1_path_out[-1]] != 0:
                                    nrb_out_v_ = list(G_direct_v.out_edges(nbunch=seg1_path_out[-1]))[0][1]
                                    seg_v_ = find_seg_passlabel_new(branchnode=seg1_path_out[-1], nbr=nrb_out_v_,
                                                                    graph=g4_not_cut, close_disc=close_disc)
                                    score_v_next = 0
                                    score_a_next = 0
                                    for node in seg_v_:
                                        if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                                            score_v_next += 1
                                        if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                                            score_a_next += 1
                                    seg_v_all_path_out = find_the_path_new(seg1_path_out[-1], nrb_out_v_,
                                                                           G_direct_allseg, g4_not_cut)
                                    seg_v_path_out = find_the_path_new(seg1_path_out[-1], nrb_out_v_, G_direct,
                                                                       g4_not_cut)
                                    if score_v_next > score_a_next:
                                        for node1_ in seg_v_path_out[1: -1]:
                                            G_direct_v.remove_node(node1_)
                                        nx.add_path(G_direct_a, seg_v_path_out)
                                        for node1_ in seg_v_all_path_out[1: -1]:
                                            G_direct_v_allseg.remove_node(node1_)
                                        nx.add_path(G_direct_a_allseg, seg_v_all_path_out)
                                if G_direct.out_degree[seg1_path_out[-1]] == G_direct_a.out_degree[seg1_path_out[-1]]:
                                    if seg1_path_out[1] in G_direct_v:
                                        for node1_ in seg1_path_out[1: -1]:
                                            G_direct_v.remove_node(node1_)
                                        nx.add_path(G_direct_a, seg1_path_out)
                                        for node1_ in seg1all_path_out[1: -1]:
                                            G_direct_v_allseg.remove_node(node1_)
                                        nx.add_path(G_direct_a_allseg, seg1all_path_out)
                        else:
                            if seg1_path_out[1] in G_direct_a:
                                for node1_ in seg1_path_out[1: -1]:
                                    G_direct_a.remove_node(node1_)
                                nx.add_path(G_direct_v, seg1_path_out)
                                for node1_ in seg1all_path_out[1: -1]:
                                    G_direct_a_allseg.remove_node(node1_)
                                nx.add_path(G_direct_v_allseg, seg1all_path_out)
                            if seg2_path_out[2] in G_direct_a:
                                continue
                            else:
                                for node1_ in seg2_path_out[1: -1]:
                                    G_direct_v.remove_node(node1_)
                                nx.add_path(G_direct_a, seg2_path_out)
                                for node1_ in seg2all_path_out[1: -1]:
                                    G_direct_v_allseg.remove_node(node1_)
                                nx.add_path(G_direct_a_allseg, seg2all_path_out)

                    else:
                        seg1_path_in = find_the_path_new(segin_1[-1], segin_1[-2], G_direct, g4_not_cut)
                        seg1_path_out = find_the_path_new(node1, nrbout_1, G_direct, g4_not_cut)
                        seg1all_path_in = find_the_path_new(segin_1[-1], segin_1[-2], G_direct_allseg, g4_not_cut)
                        seg1all_path_out = find_the_path_new(node1, nrbout_1, G_direct_allseg, g4_not_cut)
                        if seg1_path_in[-2] in G_direct_a:
                            for node1_ in seg1_path_in[1: -1]:
                                G_direct_a.remove_node(node1_)
                            nx.add_path(G_direct_v, seg1_path_in)
                            for node1_ in seg1all_path_in[1: -1]:
                                G_direct_a_allseg.remove_node(node1_)
                            nx.add_path(G_direct_v_allseg, seg1all_path_in)
                        seg2_path_out = find_the_path_new(node1, nrbout_2, G_direct, g4_not_cut)
                        if seg2_path_out[1] in G_direct_v and seg2_path_out[-1] == seg1_path_out[-1]:
                            continue
                        seg2_path_in = find_the_path_new(segin_2[-1], segin_2[-2], G_direct, g4_not_cut)
                        if seg1_path_in[0] == seg2_path_in[0] or seg1_path_out[-1] == seg2_path_out[-1]:
                            continue

                        seg2all_path_in = find_the_path_new(segin_2[-1], segin_2[-2], G_direct_allseg, g4_not_cut)
                        seg2all_path_out = find_the_path_new(node1, nrbout_2, G_direct_allseg, g4_not_cut)
                        seg2all_path_in_ = [seg2all_path_in[0], seg2all_path_in[1],
                                            seg2all_path_in[int(len(seg2all_path_in) / 2)], seg2all_path_in[-2],
                                            seg2all_path_in[-1]]
                        seg2all_path_out_ = [seg2all_path_out[0], seg2all_path_out[1],
                                             seg2all_path_out[int(len(seg2all_path_out) / 2)], seg2all_path_out[-2],
                                             seg2all_path_out[-1]]
                        seg1all_path_in_ = [seg1all_path_in[0], seg1all_path_in[1],
                                            seg1all_path_in[int(len(seg1all_path_in) / 2)], seg1all_path_in[-2],
                                            seg1all_path_in[-1]]
                        seg1all_path_out_ = [seg1all_path_out[0], seg1all_path_out[1],
                                             seg1all_path_out[int(len(seg1all_path_out) / 2)], seg1all_path_out[-2],
                                             seg1all_path_out[-1]]
                        (anglein, nrbin) = segs_angle_othertwosege1(seg2all_path_in_, seg1all_path_in_)
                        (angleout, nrbout) = segs_angle_othertwosege(seg2all_path_out_, seg1all_path_out_)
                        if angleout > 20 or anglein > 20:
                            if seg1_path_out[2] in G_direct_a:
                                for node1_ in seg1_path_out[1: -1]:
                                    G_direct_a.remove_node(node1_)
                                nx.add_path(G_direct_v, seg1_path_out)
                                for node1_ in seg1all_path_out[1: -1]:
                                    G_direct_a_allseg.remove_node(node1_)
                                nx.add_path(G_direct_v_allseg, seg1all_path_out)
                            seg2all_path_out = find_the_path_new(node1, nrbout_2, G_direct_allseg, g4_not_cut)
                            if seg2_path_out[1] in G_direct_v:
                                a1, a2, score_v_lastINv, score_a_lastINa = find_score(bn_class, G_direct_new,
                                                                                      G_direct_v,
                                                                                      G_direct_a, seg1_path_out[0], g4,
                                                                                      seg2all_path_out,
                                                                                      seg2all_path_out[1], close_disc, image_label, Disc_size)
                                if a1 > a2:
                                    for node1_ in seg2_path_out[1: -1]:
                                        G_direct_v.remove_node(node1_)
                                    nx.add_path(G_direct_a, seg2_path_out)
                                    for node1_ in seg2all_path_out[1: -1]:
                                        G_direct_v_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_a_allseg, seg2all_path_out)
                                if G_direct_v.out_degree[seg2_path_out[-1]] != 0:
                                    nrb_out_v_ = list(G_direct_v.out_edges(nbunch=seg2_path_out[-1]))[0][1]
                                    seg_v_ = find_seg_passlabel_new(branchnode=seg2_path_out[-1], nbr=nrb_out_v_,
                                                                    graph=g4_not_cut, close_disc=close_disc)
                                    score_v_next = 0
                                    score_a_next = 0
                                    for node in seg_v_:
                                        if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                                            score_v_next += 1
                                        if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                                            score_a_next += 1
                                    seg_v_all_path_out = find_the_path_new(seg2_path_out[-1], nrb_out_v_,
                                                                           G_direct_allseg, g4_not_cut)
                                    seg_v_path_out = find_the_path_new(seg2_path_out[-1], nrb_out_a_, G_direct,
                                                                       g4_not_cut)
                                    if score_v_next > score_a_next:
                                        for node1_ in seg_v_path_out[1: -1]:
                                            G_direct_v.remove_node(node1_)
                                        nx.add_path(G_direct_a, seg_v_path_out)
                                        for node1_ in seg_v_all_path_out[1: -1]:
                                            G_direct_v_allseg.remove_node(node1_)
                                        nx.add_path(G_direct_a_allseg, seg_v_all_path_out)
                                if G_direct.out_degree[seg2_path_out[-1]] == G_direct_a.out_degree[seg2_path_out[-1]]:
                                    if seg2_path_out[1] in G_direct_v:
                                        for node1_ in seg2_path_out[1: -1]:
                                            G_direct_v.remove_node(node1_)
                                        nx.add_path(G_direct_a, seg2_path_out)
                                        for node1_ in seg2all_path_out[1: -1]:
                                            G_direct_v_allseg.remove_node(node1_)
                                        nx.add_path(G_direct_a_allseg, seg2all_path_out)
                        else:
                            if seg2_path_out[1] in G_direct_a:
                                for node1_ in seg2_path_out[1: -1]:
                                    G_direct_a.remove_node(node1_)
                                nx.add_path(G_direct_v, seg2_path_out)
                                for node1_ in seg2all_path_out[1: -1]:
                                    G_direct_a_allseg.remove_node(node1_)
                                nx.add_path(G_direct_v_allseg, seg2all_path_out)
                            if seg1_path_out[2] in G_direct_a:
                                continue
                            else:
                                for node1_ in seg1_path_out[1: -1]:
                                    G_direct_v.remove_node(node1_)
                                nx.add_path(G_direct_a, seg1_path_out)
                                for node1_ in seg1all_path_out[1: -1]:
                                    G_direct_v_allseg.remove_node(node1_)
                                nx.add_path(G_direct_a_allseg, seg1all_path_out)

        if in_degree_a == 2 or in_degree_v == 2:
            if bn_class == (3, 2) or bn_class == (3, 3, 1) or bn_class == (3, 3, 2):
                nrb_in__ = []
                for edge in list(G_direct_allseg.in_edges(nbunch=node1)):
                    nrb_in__.append(edge[0])
                seg1_ = find_seg_passlabel_new(branchnode=node1, nbr=nrb_in__[0], graph=g4_not_cut, close_disc=close_disc)
                seg2_ = find_seg_passlabel_new(branchnode=node1, nbr=nrb_in__[1], graph=g4_not_cut, close_disc=close_disc)
                if seg1_[-1] == seg2_[-1]:
                    score_segin1_v = 0
                    score_segin1_a = 0
                    score_segin2_v = 0
                    score_segin2_a = 0
                    for node in seg1_:
                        if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                            score_segin1_v += 1
                        if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                            score_segin1_a += 1
                    for node in seg2_:
                        if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                            score_segin2_v += 1
                        if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                            score_segin2_a += 1
                    p_a1 = score_segin1_a/(score_segin1_v + score_segin1_a + 1e-6)
                    p_v1 = score_segin1_v/(score_segin1_v + score_segin1_a + 1e-6)
                    p_a2 = score_segin2_a/(score_segin2_v + score_segin2_a + 1e-6)
                    p_v2 = score_segin2_v / (score_segin2_v + score_segin2_a + 1e-6)
                    if (p_a1 > 0.8 and p_a2 > 0.8) or (p_v1 > 0.8 and p_v2 > 0.8):
                        continue

                nrb_out = list(G_direct.out_edges(nbunch=node1))[0][1]
                seg_out = find_seg_passlabel_new(branchnode=node1, nbr=nrb_out, graph=g4_not_cut, close_disc=close_disc)

                out_degree = (
                    0 if len({G_direct.out_degree(nbunch=seg_out[-1])}) == 0 else G_direct.out_degree(nbunch=seg_out[-1]))
                if out_degree == 0 or out_degree == 1:
                    nrb_in = []
                    for edge in list(G_direct.in_edges(nbunch=node1)):
                        nrb_in.append(edge[0])
                    segin_1 = find_seg_passlabel_new(branchnode=node1, nbr=nrb_in[0], graph=g4_not_cut,
                                                     close_disc=close_disc)
                    segin_2 = find_seg_passlabel_new(branchnode=node1, nbr=nrb_in[1], graph=g4_not_cut,
                                                     close_disc=close_disc)
                    score_segin1_v = 0
                    score_segin1_a = 0
                    score_segin2_v = 0
                    score_segin2_a = 0
                    nrb_1_lastin = []
                    nrb_2_lastin = []
                    seg_1_lastin = []
                    seg_2_lastin = []
                    score_segin1_lastin_v = 0
                    score_segin1_lastin_a = 0
                    score_segin2_lastin_v = 0
                    score_segin2_lastin_a = 0
                    for edge in list(G_direct.in_edges(nbunch=segin_1[-1])):
                        nrb_1_lastin.append(edge[0])
                    for edge in list(G_direct.in_edges(nbunch=segin_2[-1])):
                        nrb_2_lastin.append(edge[0])
                    for edge in list(g4_not_cut.edges(segin_1[-1])):
                        seg = find_seg_passlabel_new(branchnode=segin_1[-1], nbr=edge[1], graph=g4_not_cut, close_disc=close_disc)
                        if seg[1] in nrb_1_lastin:
                            seg_1_lastin.append(seg)
                    for edge in list(g4_not_cut.edges(segin_2[-1])):
                        seg = find_seg_passlabel_new(branchnode=segin_2[-1], nbr=edge[1], graph=g4_not_cut, close_disc=close_disc)
                        if seg[1] in nrb_2_lastin:
                            seg_2_lastin.append(seg)
                    proportion_v_segin1_lastin = 0
                    proportion_a_segin1_lastin = 0
                    for i in range(0, len(seg_1_lastin)):
                        seg = seg_1_lastin[i]
                        for node in seg:
                            if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                                score_segin1_lastin_v += 1
                            if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                                score_segin1_lastin_a += 1
                        proportion_v_segin1_lastin += score_segin1_lastin_v / len(seg)
                        proportion_a_segin1_lastin += score_segin1_lastin_a / len(seg)
                    proportion_v_segin1_lastin = proportion_v_segin1_lastin/(len(seg_1_lastin)+1e-8)
                    proportion_a_segin1_lastin = proportion_a_segin1_lastin/(len(seg_1_lastin)+1e-8)
                    proportion_v_segin2_lastin = 0
                    proportion_a_segin2_lastin = 0
                    for i in range(0, len(seg_2_lastin)):
                        seg = seg_2_lastin[i]
                        for node in seg:
                            if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                                score_segin2_lastin_v += 1
                            if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                                score_segin2_lastin_a += 1
                        proportion_v_segin2_lastin += score_segin2_lastin_v / len(seg)
                        proportion_a_segin2_lastin += score_segin2_lastin_a / len(seg)
                    proportion_v_segin2_lastin = proportion_v_segin2_lastin / (len(seg_2_lastin)+1e-8)
                    proportion_a_segin2_lastin = proportion_a_segin2_lastin / (len(seg_2_lastin)+1e-8)

                    for node in segin_1:
                        if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                            score_segin1_v += 1
                        if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                            score_segin1_a += 1
                    proportion_v_segin1 = score_segin1_v / len(segin_1)
                    proportion_a_segin1 = score_segin1_a / len(segin_1)
                    for node in segin_2:
                        if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                            score_segin2_v += 1
                        if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                            score_segin2_a += 1
                    proportion_v_segin2 = score_segin2_v / len(segin_2)
                    proportion_a_segin2 = score_segin2_a / len(segin_2)
                    seg1_v = proportion_v_segin1 + proportion_v_segin1_lastin
                    seg1_a = proportion_a_segin1 + proportion_a_segin1_lastin
                    seg2_v = proportion_v_segin2 + proportion_v_segin2_lastin
                    seg2_a = proportion_a_segin2 + proportion_a_segin2_lastin
                    if in_degree_v == 2:
                        if seg1_v > seg2_v:
                            seg2_path_in = find_the_path_new(segin_2[-1], segin_2[-2], G_direct, g4_not_cut)
                            seg2all_path_in = find_the_path_new(segin_2[-1], segin_2[-2], G_direct_allseg, g4_not_cut)
                            if seg2_path_in[-2] in G_direct_v:
                                for node1_ in seg2_path_in[1: -1]:
                                    G_direct_v.remove_node(node1_)
                                nx.add_path(G_direct_a, seg2_path_in)
                                for node1_ in seg2all_path_in[1: -1]:
                                    G_direct_v_allseg.remove_node(node1_)
                                nx.add_path(G_direct_a_allseg, seg2all_path_in)
                        else:
                            seg1_path_in = find_the_path_new(segin_1[-1], segin_1[-2], G_direct, g4_not_cut)
                            seg1all_path_in = find_the_path_new(segin_1[-1], segin_1[-2], G_direct_allseg, g4_not_cut)
                            if seg1_path_in[-2] in G_direct_v:
                                for node1_ in seg1_path_in[1: -1]:
                                    G_direct_v.remove_node(node1_)
                                nx.add_path(G_direct_a, seg1_path_in)
                                for node1_ in seg1all_path_in[1: -1]:
                                    G_direct_v_allseg.remove_node(node1_)
                                nx.add_path(G_direct_a_allseg, seg1all_path_in)
                    if in_degree_a == 2:
                        if seg1_a > seg2_a:
                            seg2_path_in = find_the_path_new(segin_2[-1], segin_2[-2], G_direct, g4_not_cut)
                            seg2all_path_in = find_the_path_new(segin_2[-1], segin_2[-2], G_direct_allseg, g4_not_cut)
                            if seg2_path_in[-2] in G_direct_a:
                                for node1_ in seg2_path_in[1: -1]:
                                    G_direct_a.remove_node(node1_)
                                nx.add_path(G_direct_v, seg2_path_in)
                                for node1_ in seg2all_path_in[1: -1]:
                                    G_direct_a_allseg.remove_node(node1_)
                                nx.add_path(G_direct_v_allseg, seg2all_path_in)
                        else:
                            seg1_path_in = find_the_path_new(segin_1[-1], segin_1[-2], G_direct, g4_not_cut)
                            seg1all_path_in = find_the_path_new(segin_1[-1], segin_1[-2], G_direct_allseg, g4_not_cut)
                            if seg1_path_in[-2] in G_direct_a:
                                for node1_ in seg1_path_in[1: -1]:
                                    G_direct_a.remove_node(node1_)
                                nx.add_path(G_direct_v, seg1_path_in)
                                for node1_ in seg1all_path_in[1: -1]:
                                    G_direct_a_allseg.remove_node(node1_)
                                nx.add_path(G_direct_v_allseg, seg1all_path_in)
                    segin_1_ = [segin_1[0], segin_1[1], segin_1[int(len(segin_1) / 2)], segin_1[-2], segin_1[-1]]
                    segin_2_ = [segin_2[0], segin_2[1], segin_2[int(len(segin_2) / 2)], segin_2[-2], segin_2[-1]]
                    seg_out = [seg_out[0], seg_out[1], seg_out[int(len(seg_out) / 2)], seg_out[-2], seg_out[-1]]
                    (angle1, nrb) = segs_angle_othertwosege(segin_1_, seg_out)
                    (angle2, nrb) = segs_angle_othertwosege(segin_2_, seg_out)
                    if angle1 > angle2:
                        if segin_1[1] in G_direct_a_allseg and seg_out[1] not in G_direct_a_allseg:
                            seg_path_out = find_the_path_new(seg_out[0], nrb_out, G_direct, g4_not_cut)
                            segall_path_out = find_the_path_new(seg_out[0], nrb_out, G_direct_allseg, g4_not_cut)
                            for node1_ in seg_path_out[1: -1]:
                                G_direct_v.remove_node(node1_)
                            nx.add_path(G_direct_a, seg_path_out)
                            for node1_ in segall_path_out[1: -1]:
                                G_direct_v_allseg.remove_node(node1_)
                            nx.add_path(G_direct_a_allseg, segall_path_out)
                        if segin_1[1] in G_direct_v_allseg and seg_out[1] not in G_direct_v_allseg:
                            seg_path_out = find_the_path_new(seg_out[0], nrb_out, G_direct, g4_not_cut)
                            segall_path_out = find_the_path_new(seg_out[0], nrb_out, G_direct_allseg, g4_not_cut)
                            for node1_ in seg_path_out[1: -1]:
                                G_direct_a.remove_node(node1_)
                            nx.add_path(G_direct_v, seg_path_out)
                            for node1_ in segall_path_out[1: -1]:
                                G_direct_a_allseg.remove_node(node1_)
                            nx.add_path(G_direct_v_allseg, segall_path_out)
                    else:
                        if segin_2[1] in G_direct_a_allseg and seg_out[1] not in G_direct_a_allseg:
                            seg_path_out = find_the_path_new(seg_out[0], nrb_out, G_direct, g4_not_cut)
                            segall_path_out = find_the_path_new(seg_out[0], nrb_out, G_direct_allseg, g4_not_cut)
                            for node1_ in seg_path_out[1: -1]:
                                G_direct_v.remove_node(node1_)
                            nx.add_path(G_direct_a, seg_path_out)
                            for node1_ in segall_path_out[1: -1]:
                                G_direct_v_allseg.remove_node(node1_)
                            nx.add_path(G_direct_a_allseg, segall_path_out)
                        if segin_2[1] in G_direct_v_allseg and seg_out[1] not in G_direct_v_allseg:
                            seg_path_out = find_the_path_new(seg_out[0], nrb_out, G_direct, g4_not_cut)
                            segall_path_out = find_the_path_new(seg_out[0], nrb_out, G_direct_allseg, g4_not_cut)
                            for node1_ in seg_path_out[1: -1]:
                                G_direct_a.remove_node(node1_)
                            nx.add_path(G_direct_v, seg_path_out)
                            for node1_ in segall_path_out[1: -1]:
                                G_direct_a_allseg.remove_node(node1_)
                            nx.add_path(G_direct_v_allseg, segall_path_out)

                else:
                    [(nrbin_1, nrbout_1), (nrbin_2, nrbout_2), (nrb_1, nrb_2)] = find_twoClasses_cp_new1(G_direct,
                                                                                                         node1, g4,
                                                                                                         close_disc)
                    segin_1 = find_seg_passlabel(branchnode=node1, nbr=nrbin_1, graph=g4, close_disc=close_disc)
                    segout_1 = find_seg_passlabel(branchnode=seg_out[-1], nbr=nrbout_1, graph=g4, close_disc=close_disc)
                    segin_2 = find_seg_passlabel(branchnode=node1, nbr=nrbin_2, graph=g4, close_disc=close_disc)
                    segout_2 = find_seg_passlabel(branchnode=seg_out[-1], nbr=nrbout_2, graph=g4, close_disc=close_disc)
                    score_segin1_v = 0
                    score_segin1_a = 0
                    score_segout1_v = 0
                    score_segout1_a = 0
                    score_segin2_v = 0
                    score_segin2_a = 0
                    score_segout2_v = 0
                    score_segout2_a = 0
                    for node in segin_1:
                        if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                            score_segin1_v += 1
                        if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                            score_segin1_a += 1
                    proportion_v_segin1 = score_segin1_v / (score_segin1_v+score_segin1_a + 1e-6)
                    proportion_a_segin1 = score_segin1_a / (score_segin1_v+score_segin1_a + 1e-6)
                    for node in segout_1:
                        if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                            score_segout1_v += 1
                        if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                            score_segout1_a += 1
                    proportion_v_segout1 = score_segout1_v / (score_segout1_v+score_segout1_a + 1e-6)
                    proportion_a_segout1 = score_segout1_a / (score_segout1_v+score_segout1_a + 1e-6)
                    for node in segin_2:
                        if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                            score_segin2_v += 1
                        if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                            score_segin2_a += 1
                    proportion_v_segin2 = score_segin2_v / (score_segin2_v+score_segin2_a + 1e-6)
                    proportion_a_segin2 = score_segin2_a / (score_segin2_v+score_segin2_a + 1e-6)
                    for node in segout_2:
                        if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                            score_segout2_v += 1
                        if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                            score_segout2_a += 1
                    proportion_v_segout2 = score_segout2_v / (score_segout2_v+score_segout2_a + 1e-6)
                    proportion_a_segout2 = score_segout2_a / (score_segout2_v+score_segout2_a + 1e-6)
                    abs_dist = math.sqrt((seg_out[-1][0] - seg_out[0][0]) ** 2 + (seg_out[-1][1] - seg_out[0][1]) ** 2)
                    if (len(seg_out) > 12 * 2.5 or abs_dist > 15 * 2.5 or segout_1[-1] == segout_2[-1]) and segin_1[-1] != segin_2[-1]:
                        seg1_v = proportion_v_segin1
                        seg1_a = proportion_a_segin1
                        seg2_v = proportion_v_segin2
                        seg2_a = proportion_a_segin2
                        if in_degree_v == 2:
                            if seg1_v == seg2_v:
                                nrbin3_segin1 = []
                                nrbin3_segin2 = []
                                for edge in list(G_direct.in_edges(nbunch=segin_1[-1])):
                                    nrbin3_segin1.append(edge[0])
                                angle1 = 0
                                for node_nbr in nrbin3_segin1:
                                    seg = find_seg_passlabel(branchnode=segin_1[-1], nbr=node_nbr, graph=g4,
                                                             close_disc=close_disc)
                                    (angle, nrb0) = segs_angle(seg, segin_1)
                                    if angle > angle1:
                                        [angle1, nrb_new] = [angle, nbr0]
                                seg_in1_last = find_seg_passlabel(branchnode=segin_1[-1], nbr=nrb_new, graph=g4,
                                                                  close_disc=close_disc)
                                for edge in list(G_direct.in_edges(nbunch=segin_2[-1])):
                                    nrbin3_segin2.append(edge[0])
                                angle1 = 0
                                for node_nbr in nrbin3_segin2:
                                    seg = find_seg_passlabel(branchnode=segin_2[-1], nbr=node_nbr, graph=g4,
                                                             close_disc=close_disc)
                                    (angle, nrb0) = segs_angle(seg, segin_2)
                                    if angle > angle1:
                                        [angle1, nrb_new] = [angle, nbr0]
                                seg_in2_last = find_seg_passlabel(branchnode=segin_2[-1], nbr=nrb_new, graph=g4,
                                                                  close_disc=close_disc)
                                score_seg_in1_last_v = 0
                                score_seg_in1_last_a = 0
                                score_seg_in2_last_v = 0
                                score_seg_in2_last_a = 0
                                for node in seg_in1_last:
                                    if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                                        score_seg_in1_last_v += 1
                                    if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                                        score_seg_in1_last_a += 1
                                for node in seg_in2_last:
                                    if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                                        score_seg_in2_last_v += 1
                                    if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                                        score_seg_in2_last_a += 1
                                proportion_v_seg_in1_last = score_seg_in1_last_v / (score_seg_in1_last_v + score_seg_in1_last_a + 1e-6)
                                proportion_a_seg_in1_last = score_seg_in1_last_a / (score_seg_in1_last_v + score_seg_in1_last_a + 1e-6)
                                proportion_v_seg_in2_last = score_seg_in2_last_v / (score_seg_in2_last_v + score_seg_in2_last_a + 1e-6)
                                proportion_a_seg_in2_last = score_seg_in2_last_a / (score_seg_in2_last_v + score_seg_in2_last_a + 1e-6)
                                seg1_v = seg1_v + proportion_v_seg_in1_last
                                seg1_a = seg1_a + proportion_a_seg_in1_last
                                seg2_v = seg2_v + proportion_v_seg_in2_last
                                seg2_a = seg2_a + proportion_a_seg_in2_last
                        if in_degree_v == 2:
                            if seg1_v > seg2_v:
                                seg2_path_in = find_the_path(segin_2[-1], segin_2[-2], G_direct)
                                seg2all_path_in = find_the_path(segin_2[-1], segin_2[-2], G_direct_allseg)
                                if seg2_path_in[-2] in G_direct_v:
                                    for node1_ in seg2_path_in[1: -1]:
                                        G_direct_v.remove_node(node1_)
                                    nx.add_path(G_direct_a, seg2_path_in)
                                    for node1_ in seg2all_path_in[1: -1]:
                                        G_direct_v_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_a_allseg, seg2all_path_in)
                            else:
                                seg1_path_in = find_the_path(segin_1[-1], segin_1[-2], G_direct)
                                seg1all_path_in = find_the_path(segin_1[-1], segin_1[-2], G_direct_allseg)
                                if seg1_path_in[-2] in G_direct_v:
                                    for node1_ in seg1_path_in[1: -1]:
                                        G_direct_v.remove_node(node1_)
                                    nx.add_path(G_direct_a, seg1_path_in)
                                    for node1_ in seg1all_path_in[1: -1]:
                                        G_direct_v_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_a_allseg, seg1all_path_in)
                        if in_degree_a == 2:
                            if seg1_a > seg2_a:
                                seg2_path_in = find_the_path(segin_2[-1], segin_2[-2], G_direct)
                                seg2all_path_in = find_the_path(segin_2[-1], segin_2[-2], G_direct_allseg)
                                if seg2_path_in[-2] in G_direct_a:
                                    for node1_ in seg2_path_in[1: -1]:
                                        G_direct_a.remove_node(node1_)
                                    nx.add_path(G_direct_v, seg2_path_in)
                                    for node1_ in seg2all_path_in[1: -1]:
                                        G_direct_a_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_v_allseg, seg2all_path_in)
                            else:
                                seg1_path_in = find_the_path(segin_1[-1], segin_1[-2], G_direct)
                                seg1all_path_in = find_the_path(segin_1[-1], segin_1[-2], G_direct_allseg)
                                if seg1_path_in[-2] in G_direct_a:
                                    for node1_ in seg1_path_in[1: -1]:
                                        G_direct_a.remove_node(node1_)
                                    nx.add_path(G_direct_v, seg1_path_in)
                                    for node1_ in seg1all_path_in[1: -1]:
                                        G_direct_a_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_v_allseg, seg1all_path_in)
                    else:
                        seg1_v = proportion_v_segin1 + proportion_v_segout1
                        seg1_a = proportion_a_segin1 + proportion_a_segout1
                        seg2_v = proportion_v_segin2 + proportion_v_segout2
                        seg2_a = proportion_a_segin2 + proportion_a_segout2
                        if in_degree_v == 2:
                            if seg1_v > seg2_v:
                                seg2_path_in = find_the_path(segin_2[-1], segin_2[-2], G_direct)
                                seg2_path_out = find_the_path(seg_out[-1], nrbout_2, G_direct)
                                seg2all_path_in = find_the_path(segin_2[-1], segin_2[-2], G_direct_allseg)
                                seg2all_path_out = find_the_path(seg_out[-1], nrbout_2, G_direct_allseg)
                                if seg2_path_in[-2] in G_direct_v:
                                    for node1_ in seg2_path_in[1: -1]:
                                        G_direct_v.remove_node(node1_)
                                    nx.add_path(G_direct_a, seg2_path_in)
                                    for node1_ in seg2all_path_in[1: -1]:
                                        G_direct_v_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_a_allseg, seg2all_path_in)
                                if seg2_path_out[2] in G_direct_v:
                                    for node1_ in seg2_path_out[1: -1]:
                                        G_direct_v.remove_node(node1_)
                                    nx.add_path(G_direct_a, seg2_path_out)
                                    for node1_ in seg2all_path_out[1: -1]:
                                        G_direct_v_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_a_allseg, seg2all_path_out)
                            else:
                                seg1_path_in = find_the_path(segin_1[-1], segin_1[-2], G_direct)
                                seg1_path_out = find_the_path(seg_out[-1], nrbout_1, G_direct)
                                seg1all_path_in = find_the_path(segin_1[-1], segin_1[-2], G_direct_allseg)
                                seg1all_path_out = find_the_path(seg_out[-1], nrbout_1, G_direct_allseg)
                                if seg1_path_in[-2] in G_direct_v:
                                    for node1_ in seg1_path_in[1: -1]:
                                        G_direct_v.remove_node(node1_)
                                    nx.add_path(G_direct_a, seg1_path_in)
                                    for node1_ in seg1all_path_in[1: -1]:
                                        G_direct_v_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_a_allseg, seg1all_path_in)
                                if seg1_path_out[2] in G_direct_v:
                                    for node1_ in seg1_path_out[1: -1]:
                                        G_direct_v.remove_node(node1_)
                                    nx.add_path(G_direct_a, seg1_path_out)
                                    for node1_ in seg1all_path_out[1: -1]:
                                        G_direct_v_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_a_allseg, seg1all_path_out)
                        if in_degree_a == 2:
                            if seg1_a > seg2_a:
                                seg2_path_in = find_the_path(segin_2[-1], segin_2[-2], G_direct)
                                seg2_path_out = find_the_path(seg_out[-1], nrbout_2, G_direct)
                                seg2all_path_in = find_the_path(segin_2[-1], segin_2[-2], G_direct_allseg)
                                seg2all_path_out = find_the_path(seg_out[-1], nrbout_2, G_direct_allseg)
                                if seg2_path_in[-2] in G_direct_a:
                                    for node1_ in seg2_path_in[1: -1]:
                                        G_direct_a.remove_node(node1_)
                                    nx.add_path(G_direct_v, seg2_path_in)
                                    for node1_ in seg2all_path_in[1: -1]:
                                        G_direct_a_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_v_allseg, seg2all_path_in)
                                if seg2_path_out[2] in G_direct_a:
                                    for node1_ in seg2_path_out[1: -1]:
                                        G_direct_a.remove_node(node1_)
                                    nx.add_path(G_direct_v, seg2_path_out)
                                    for node1_ in seg2all_path_out[1: -1]:
                                        G_direct_a_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_v_allseg, seg2all_path_out)
                            else:
                                seg1_path_in = find_the_path(segin_1[-1], segin_1[-2], G_direct)
                                seg1_path_out = find_the_path(seg_out[-1], nrbout_1, G_direct)
                                seg1all_path_in = find_the_path(segin_1[-1], segin_1[-2], G_direct_allseg)
                                seg1all_path_out = find_the_path(seg_out[-1], nrbout_1, G_direct_allseg)
                                if seg1_path_in[-2] in G_direct_a:
                                    for node1_ in seg1_path_in[1: -1]:
                                        G_direct_a.remove_node(node1_)
                                    nx.add_path(G_direct_v, seg1_path_in)
                                    for node1_ in seg1all_path_in[1: -1]:
                                        G_direct_a_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_v_allseg, seg1all_path_in)
                                if seg1_path_out[2] in G_direct_a:
                                    for node1_ in seg1_path_out[1: -1]:
                                        G_direct_a.remove_node(node1_)
                                    nx.add_path(G_direct_v, seg1_path_out)
                                    for node1_ in seg1all_path_out[1: -1]:
                                        G_direct_a_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_v_allseg, seg1all_path_out)

        in_degree_a_1 = (0 if len({G_direct_a.in_degree(nbunch=node1)}) == 0 else G_direct_a.in_degree(nbunch=node1))
        in_degree_v_2 = (0 if len({G_direct_v.in_degree(nbunch=node1)}) == 0 else G_direct_v.in_degree(nbunch=node1))
        if (in_degree_a == 2 and in_degree_a_1 == 1) or (in_degree_v == 2 and in_degree_v_2 == 1):
            nrb_in = []
            for edge in list(G_direct.in_edges(nbunch=node1)):
                nrb_in.append(edge[0])
            segin_1 = find_seg_passlabel_new(branchnode=node1, nbr=nrb_in[0], graph=g4_not_cut, close_disc=close_disc)
            segin_2 = find_seg_passlabel_new(branchnode=node1, nbr=nrb_in[1], graph=g4_not_cut, close_disc=close_disc)
            segin_1_path = find_the_path_new(segin_1[-1], segin_1[-2], G_direct, g4_not_cut)
            if in_degree_a == 2:
                if segin_1_path[1] not in G_direct_a:
                    segin_correct = segin_1
                else:
                    segin_correct = segin_2
            else:
                if segin_1_path[1] not in G_direct_v:
                    segin_correct = segin_1
                else:
                    segin_correct = segin_2
            while 1:
                i, j = 0, 0
                seg_next = segin_correct
                in_degree_start = (0 if len({G_direct_new.in_degree(nbunch=seg_next[-1])}) == 0 else G_direct_new.in_degree(nbunch=seg_next[-1]))
                dist_nextseg_start = math.sqrt((seg_next[-1][0] - H[0]) ** 2 + (seg_next[-1][1] - H[1]) ** 2)
                out_degreeV_star = (0 if len({G_direct_v.out_degree(nbunch=seg_next[-1])}) == 0 else G_direct_v.out_degree(nbunch=seg_next[-1]))
                out_degreeA_star = (0 if len({G_direct_a.out_degree(nbunch=seg_next[-1])}) == 0 else G_direct_a.out_degree(nbunch=seg_next[-1]))

                if in_degree_start >= 2 or in_degree_start == 0 or dist_nextseg_start < Disc_size * 2:
                    break
                else:
                    if (seg_next[-2] in G_direct_a_allseg and out_degreeV_star == 0) or (seg_next[-2] in G_direct_v_allseg and out_degreeA_star == 0):
                        nrb_out_seg_next0 = []
                        for edge in list(G_direct_new.out_edges(nbunch=seg_next[-1])):
                            nrb_out_seg_next0.append(edge[1])
                        for node_ in nrb_out_seg_next0:
                            seg_ = find_seg_passlabel_new(branchnode=seg_next[-1], nbr=node_, graph=g4_not_cut,
                                                         close_disc=close_disc)
                            out_degreeV_ = (0 if len({G_direct_v.out_degree(nbunch=seg_[-1])}) == 0 else G_direct_v.out_degree(nbunch=seg_next[-1]))
                            out_degreeA_ = (0 if len({G_direct_a.out_degree(nbunch=seg_[-1])}) == 0 else G_direct_a.out_degree(nbunch=seg_next[-1]))
                            if seg_next[-2] not in seg_ and out_degreeV_ != 0 and out_degreeA_ != 0 and G_direct_new.in_degree[seg_[-1]] == 1:
                                break
                    seg_last_in = find_seg_passlabel_new(branchnode=seg_next[-1], nbr=list(G_direct.in_edges(nbunch=seg_next[-1]))[0][0], graph=g4_not_cut,
                                                     close_disc=close_disc)
                    bn_class_in = find_branchnode_class(G_direct_new, seg_next[-1], g4_not_cut, close_disc, H, Disc_size)
                    b1, b2, score_v_lastINv, score_a_lastINa = find_score(bn_class_in, G_direct_new, G_direct_v,
                                                                          G_direct_a, seg_next[-1], g4, seg_next,
                                                                          seg_next[-2], close_disc, image_label, H,
                                                                          Disc_size)
                    nrb_out_start = []
                    for edge in list(G_direct.out_edges(nbunch=seg_next[-1])):
                        if edge[1] not in seg_next:
                            nrb_out_start.append(edge[1])
                    score_v_next = 0
                    score_a_next = 0
                    for node in seg_next:
                        if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                            score_v_next += 1
                        if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                            score_a_next += 1
                    proportion_v_next_0 = score_v_next / (score_v_next + score_a_next + 1e-6)
                    proportion_a_next_0 = score_a_next / (score_v_next + score_a_next + 1e-6)
                    if seg_next[-2] in G_direct_a_allseg and score_a_lastINa == 0:
                        score_a_lastINa_ = score_v_lastINv
                        score_v_lastINv_ = 0
                        proportion_v_next = b1 * score_v_lastINv_ + b2 * proportion_v_next_0
                        proportion_a_next = b1 * score_a_lastINa_ + b2 * proportion_a_next_0
                        if proportion_a_next > proportion_v_next:
                            seg_last_in_path = find_the_path(seg_last_in[-1], seg_last_in[-2], G_direct)
                            segall_last_in_path = seg_last_in[::-1]
                            if seg_last_in_path[-2] in G_direct_v:
                                i = 1
                                for node1_ in seg_last_in_path[1: -1]:
                                    G_direct_v.remove_node(node1_)
                                nx.add_path(G_direct_a, seg_last_in_path)
                                for node1_ in segall_last_in_path[1: -1]:
                                    G_direct_v_allseg.remove_node(node1_)
                                nx.add_path(G_direct_a_allseg, segall_last_in_path)
                    if seg_next[-2] in G_direct_v_allseg and score_v_lastINv == 0:
                        score_v_lastINv_ = score_a_lastINa
                        score_a_lastINa_ = 0
                        proportion_v_next = b1 * score_v_lastINv_ + b2 * proportion_v_next_0
                        proportion_a_next = b1 * score_a_lastINa_ + b2 * proportion_a_next_0
                        if proportion_v_next > proportion_a_next:
                            seg_last_in_path = find_the_path(seg_last_in[-1], seg_last_in[-2], G_direct)
                            segall_last_in_path = seg_last_in[::-1]
                            if seg_last_in_path[-2] in G_direct_a:
                                j = 1
                                for node1_ in seg_last_in_path[1: -1]:
                                    G_direct_a.remove_node(node1_)
                                nx.add_path(G_direct_v, seg_last_in_path)
                                for node1_ in segall_last_in_path[1: -1]:
                                    G_direct_a_allseg.remove_node(node1_)
                                nx.add_path(G_direct_v_allseg, segall_last_in_path)
                nrb_in1 = []
                for edge in list(G_direct_new.in_edges(nbunch=seg_last_in[-1])):
                    nrb_in1.append(edge[0])
                ina = 100
                inv = 100
                num = 50
                if len(nrb_in1) == 1:
                    seg_inlast = find_seg_passlabel_new(branchnode=seg_last_in[-1], nbr=nrb_in1[0],
                                                        graph=g4_not_cut, close_disc=close_disc)
                    ina = (0 if len({G_direct_a.in_degree(nbunch=seg_inlast[-1])}) == 0 else G_direct_a.in_degree(nbunch=seg_inlast[-1]))
                    inv = (0 if len({G_direct_v.in_degree(nbunch=seg_inlast[-1])}) == 0 else G_direct_v.in_degree(nbunch=seg_inlast[-1]))
                    num = G_direct_new.in_degree[seg_inlast[-1]]
                if (score_a_lastINa == 0 and inv == num) or (
                        score_v_lastINv == 0 and ina == num) or i != 0 or j != 0:
                    for nrb in nrb_out_start:
                        seg_branch = find_seg_passlabel_new(branchnode=seg_next[-1], nbr=nrb, graph=g4_not_cut, close_disc=close_disc)
                        seg_branch_path = find_the_path(seg_branch[0], seg_branch[1], G_direct)
                        score_v_next = 0
                        score_a_next = 0
                        for node in seg_branch:
                            if (image_label[node[0], node[1], :] == [255, 0, 0]).all():
                                score_v_next += 1
                            if (image_label[node[0], node[1], :] == [0, 0, 255]).all():
                                score_a_next += 1
                        proportion_v_next_0 = score_v_next / (score_v_next + score_a_next + 1e-6)
                        proportion_a_next_0 = score_a_next / (score_v_next + score_a_next + 1e-6)
                        if i == 1:
                            score_a_lastINa_ = score_v_lastINv
                            score_v_lastINv_ = 0
                            proportion_v_next = b1 * score_v_lastINv_ + b2 * proportion_v_next_0
                            proportion_a_next = b1 * score_a_lastINa_ + b2 * proportion_a_next_0
                            if proportion_a_next > proportion_v_next:
                                if seg_branch_path[-2] in G_direct_v:
                                    i = 1
                                    for node1_ in seg_branch_path[1: -1]:
                                        G_direct_v.remove_node(node1_)
                                    nx.add_path(G_direct_a, seg_branch_path)
                                    for node1_ in seg_branch[1: -1]:
                                        G_direct_v_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_a_allseg, seg_branch)
                        elif j == 1:
                            score_v_lastINv_ = score_a_lastINa
                            score_a_lastINa_ = 0
                            proportion_v_next = b1 * score_v_lastINv_ + b2 * proportion_v_next_0
                            proportion_a_next = b1 * score_a_lastINa_ + b2 * proportion_a_next_0
                            if proportion_v_next > proportion_a_next:
                                if seg_branch_path[-2] in G_direct_a:
                                    j = 1
                                    for node1_ in seg_branch_path[1: -1]:
                                        G_direct_a.remove_node(node1_)
                                    nx.add_path(G_direct_v, seg_branch_path)
                                    for node1_ in seg_branch[1: -1]:
                                        G_direct_a_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_v_allseg, seg_branch)
                        else:
                            if score_a_lastINa == 0:
                                if seg_branch_path[-2] in G_direct_a:
                                    for node1_ in seg_branch_path[1: -1]:
                                        G_direct_a.remove_node(node1_)
                                    nx.add_path(G_direct_v, seg_branch_path)
                                    for node1_ in seg_branch[1: -1]:
                                        G_direct_a_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_v_allseg, seg_branch)
                            if score_v_lastINv == 0:
                                if seg_branch_path[-2] in G_direct_v:
                                    for node1_ in seg_branch_path[1: -1]:
                                        G_direct_v.remove_node(node1_)
                                    nx.add_path(G_direct_a, seg_branch_path)
                                    for node1_ in seg_branch[1: -1]:
                                        G_direct_v_allseg.remove_node(node1_)
                                    nx.add_path(G_direct_a_allseg, seg_branch)
                segin_correct = seg_last_in

    return G_direct_v, G_direct_a, G_direct, G_direct_v_allseg, G_direct_a_allseg, G_direct_allseg
