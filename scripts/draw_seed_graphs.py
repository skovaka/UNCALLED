import sys, cv2
import numpy as np

WIN_WIDTH = 2000
WIN_HEIGHT = 1200

NODE_RADIUS = 5

node_ids = dict()

def parse_graph(infile, print_graph):
    global node_ids

    graph = dict()

    infile.readline()

    for line in infile:
        if line.strip() == "== end ==":
            break
        elif line[0] == "=":
            continue
        
        tabs = line.strip().split()
        try:
            event, length, pcount, ntype = map(int, tabs[1:5])
        except e:
            print (line.strip())
            return None
        parent = tabs[0]

        seed_node = (event == length - 1)


        if not parent in node_ids:
            node_ids[parent] = len(node_ids)

        pid = node_ids[parent]

        graph[pid] = (event, list(), seed_node)

        children = list() if len(tabs) == 5 else tabs[5:]

        for c in children:
            if not c in node_ids:
                #print (c, len(node_ids))
                node_ids[c] = len(node_ids)

            cid = node_ids[c]
            graph[pid][1].append(cid)

            if print_graph:
                print ("\t%s -> %s;" % (node_ids[parent], node_ids[c]))

    return graph

def get_win_coords(e, n, event_sizes, graph_width = 32):
    x = int( e * float(WIN_WIDTH / graph_width) )
    y = int( n * float(WIN_HEIGHT / event_sizes[e]))
    return x, y

def draw_graph(graph, img, graph_width = 32):
    
    img.fill(255)


    event_sizes = list()
    node_coords = dict()

    for n in sorted(graph):
        e, children, seed_node = graph[n]

        while e >= len(event_sizes):
            event_sizes.append(0)

        if not n in node_coords:
            node_coords[n] = (e, event_sizes[e])
            event_sizes[e] += 1

        for c in children:
            while e + 1 >= len(event_sizes):
                event_sizes.append(0)

            if not c in node_coords:
                node_coords[c] = (e+1, event_sizes[e+1])
                event_sizes[e+1] += 1

    for n in sorted(graph):
        e, children, seed_node = graph[n]

        if not seed_node:
            color = (0, 0, 0)
            size = 1

            x, y = get_win_coords(e, node_coords[n][1], event_sizes)
            cv2.circle(img, (x, y), 2, color, -1)
            
            for c in children:

                x2, y2 = get_win_coords(e+1, node_coords[c][1], event_sizes)
                cv2.line(img, (x, y), (x2, y2), color, size, cv2.LINE_AA)


    for n in sorted(graph):
        e, children, seed_node = graph[n]

        if seed_node:
            color = (0, 0, 255)
            size = 2

            x, y = get_win_coords(e, node_coords[n][1], event_sizes)
            cv2.circle(img, (x, y), 3, color, -1)

            for c in children:
                x2, y2 = get_win_coords(e+1, node_coords[c][1], event_sizes)
                cv2.line(img, (x, y), (x2, y2), color, size, cv2.LINE_AA)



    return img

if __name__ == "__main__":
    infile = open(sys.argv[1])

    if len(sys.argv) > 2:
        out_prefix = sys.argv[2]
    else:
        out_prefix = None

    img = np.ones((WIN_HEIGHT, WIN_WIDTH, 3), np.uint8)

    count = 0

    while True:
        print("Graph %d" % count)
        graph = parse_graph(infile, False)

        if len(graph) == 0:
            break
        
        draw_graph(graph, img)
        if out_prefix:
            cv2.imwrite("%s%03d.jpg" % (out_prefix, count), img)

        count += 1


