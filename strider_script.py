peds.mapper.path_reads()
from collections import Counter
def strider_walk_out(nid):
    #p holds the immediate previous node to a node in the path
    #v holds the number of times a nod appears in paths
    p={}
    v=Counter()
    if nid<0: rids=[-x for x in peds.mapper.paths_in_node[-nid]]
    else: rids=peds.mapper.paths_in_node[nid]
    for rid in rids:
        #pfw=tuple(peds.mapper.path_fw(rid,nid,True))
        #if pfw: v.update(pfw)
        v.update([x for x in peds.mapper.path_fw(rid,nid) if x!=0])
    #print(v)
    path=[nid]
    nv=ws.sdg.get_nodeview(nid)
    while True:
        last_v=0
        winner=0
        for l in nv.next():
            ln=l.node().node_id()
            if v[ln]> 5*last_v:# and v[ln]> 2:
                last_v=v[ln]
                winner=ln
            elif v[ln] >=last_v/5:
                winner = 0
                last_v = max(last_v,v[ln])
        if winner!=0 and winner not in path:
            path.append(winner)
            nv=ws.sdg.get_nodeview(winner)
        else: break
    return path

if __name__ == "__main__":
    asize=0
    anchors=set()
    for nv in ws.sdg.get_all_nodeviews():
        if nv.size()>=90 and nv.kci()>.5 and nv.kci()<1.5:
            anchors.add(nv.node_id())
            asize+=nv.size()
    print(len(anchors),'anchors found, totaling',asize,"bp")

    fw={}
    bw={}
    linear_anchors=set()
    for nid in anchors:
        fw[nid]=[strider_walk_out(nid)[1:]]
        bw[nid]=[strider_walk_out(-nid)[1:]]
        if fw[nid][0] and bw[nid][0]: linear_anchors.add(nid)
        
    print("%d / %d anchors are linear" % (len(linear_anchors),len(anchors)))

    ge=SDG.GraphEditor(ws)
    used_nodes=[]
    paths=0
    nopaths=0
    for fnid in list(linear_anchors):
        for nid in [fnid,-fnid]:
            #print("\nNode",nid)
            if nid>0: fp=fw[nid][0]
            else: fp=bw[-nid][0]
            #print("FW",fp)
            p=[]
            for x in fp:
                if abs(x) in linear_anchors:
                    p=fp[:fp.index(x)+1]
                    #print(p)
                    if x<0: op= [-x for x in fw[-x][0][::-1]]
                    else: op= [-x for x in bw[x][0][::-1]]
                    #print(op)
                    #print(p[:-1],"==",op[op.index(nid):],"??")
                    if nid not in op or p[:-1]!=op[op.index(nid)+1:]: p=[]
                    break
            #print(p)
            if p and abs(p[-1])>abs(nid):
                #print("expansion to queue",[nid]+p)
                ge.queue_path_detachment([nid]+p,True)
                try:
                    SDG.SequenceDistanceGraphPath(ws.sdg,[nid]+p).sequence()
                    paths+=1
              except:
                    nopaths+=1
                for r in p[:-1]: used_nodes.append(abs(r))
    print("%d paths and %d nopaths"%(paths,nopaths))
    ge.apply_all()

    #delete "used up" tips
    for r in range(10):
        ge=SDG.GraphEditor(ws)
        for x in used_nodes:
            try:
                nv=ws.sdg.get_nodeview(x)
                if len(nv.prev())==0 or len(nv.next())==0:
                    ge.queue_node_deletion(x)
            except: pass

        ge.apply_all()

    ws.sdg.join_all_unitigs()

