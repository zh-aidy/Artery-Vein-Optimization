import networkx as nx

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
            if image[i - 1 + hang, j - 1 + lie] == 255:
                count += 1
                if i - 1 + hang < i:
                    continue
                elif i - 1 + hang == i and j - 1 + lie < j:
                    continue
                else:
                    g.add_edge((i, j), (i - 1 + hang, j - 1 + lie))
    return count, g.edges

def mul_ske_tograph(image_zeros, h, w):
    g = nx.Graph()
    B_Gu = []
    D_g = []
    for i in range(0, h):
        for j in range(0, w):
            if image_zeros[i, j] != 0:
                g.add_node((i, j))
    for node in g.nodes:
        (i, j) = node
        k1, k2 = 0, 3
        count, Edge = neighborhood(i=i, j=j, image=image_zeros, h=h, g=g, k1=k1, k2=k2)
    for node in g.nodes:
        if g.degree[node] == 1:
            D_g.append(node)
        elif g.degree[node] > 2:
            B_Gu.append(node)
    return g, B_Gu