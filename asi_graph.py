import networkx as nt
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import itertools 
from itertools import permutations

#fonction permettant de faire un produit de deux liste (pris de itertools mais un peu modifié)
def product(*args, repeat=1):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = [tuple(pool) for pool in args] * repeat
    result = [[]]
    for pool in pools:
        result = [x+y for x in result for y in pool]
        
    for prod in result:
        yield list(prod)

#recupere l arete de départ
def get_nc(motif):
    for i in motif.edges.keys():
        if motif.edges[i]['label']!='CWW' and motif.edges[i]['label']!='B53':
            return i
    return i

#fonction permettant de comparer les labels, return la penalité de match entre deux arrêtes 
def compare_labels(a,b):
    res_range=0 #penalité de range
    if(a['long_range']!=b['long_range']):
        res_range=1
    #permet de parsec l'aret et de comparer les labels
    a=a["label"]
    b=b["label"]
    miss_match=0
    lista=a.split('+')
    listb=b.split('+')
    for i in lista:
        if len(listb)==0:
            return miss_match
        to_remove=""
        for y in listb:
            if y==i:
                to_remove=y
                break
        if to_remove=="":
            miss_match=miss_match+1
        else:
            listb.remove(to_remove)
    return miss_match,res_range

#fonction permettant de creer tout les matchs de taille 1
def get_started(g,motif,nb,nb_range):
    arr=get_nc(motif)
    #recupere l'arete de départ
    starter=[]
    for i in g.edges.keys():#pour chaque arete de la stucture
        diff,diff_range=compare_labels(g.edges[i],motif.edges[arr])#on calcule la difference entre l'arete de depart et notre arete i
        tmp_nb=0
        tmp_nb_range=0
        #attribution des penalités
        if nb_range!=-1:
            tmp_nb_range=nb_range-diff_range
        if nb!=-1:
            tmp_nb=nb-diff
        s=[(i,arr)]
        starter.append(([],s,tmp_nb,tmp_nb_range)) #on ajoute a la liste le match

    return (starter)

#fonction permettant de verifier si l'arete a deja été match avec d'autre arete de la structure pour un match donné
def check_redundancy(a,matched,just_added):
    for i in matched:
        x,y=i[0]
        if (x==a[0] and y==a[1]) or (x==a[1] and y==a[0]):
            return True
    for i in just_added:
        x,y=i[0]
        if (x==a[0] and y==a[1]) or (x==a[1] and y==a[0]):
            return True
    return False

#fonction permettant de verifier si l'arete a deja été match avec d'autre arete du motif pour un match donné
def check_redundancy_motif(a,matched,just_added):
    for i in matched:
        x,y=i[1]
        if (x==a[0] and y==a[1]) or (x==a[1] and y==a[0]):
            return True
    for i in just_added:
        x,y=i[1]
        if (x==a[0] and y==a[1]) or (x==a[1] and y==a[0]):
            return True
      
    return False

#regarde si un graph est totalement cyclique ie. si aucune arrête n'appartient pas à un cycle
def is_fully_cycled(matched):
    m=[item for sublist in [a for (a,_) in matched] for item in sublist]
    #on compte la presence de chaque noeud, si un noeud est present moins de deux fois alors le graph n est pas cyclique
    for i in set(m):
        if m.count(i)<2:
            return False
    return True

#fonction permettant de comptabiliser les aretes en plus et en moins afin de voir si la potentielle solution est une bonne ou mauvaise solution
def complete2(matched,allowed_plus,allowed_minus,sequence,motif):
    nodes_seq=[]
    nodes_mot=[]
    for a,b in matched:
        a1,a2=a
        b1,b2=b
        if a1 not in nodes_seq and b1 not in nodes_mot:
            nodes_seq.append(a1)
            nodes_mot.append(b1)
        if a2 not in nodes_seq and b2 not in nodes_mot:
            nodes_seq.append(a2)
            nodes_mot.append(b2)
    if len(nodes_mot)!=len(nodes_seq):
        return None
    nodes_seq=list(nodes_seq)
    nodes_mot=list(nodes_mot)
    for i in range(len(nodes_mot)):
        for y in range(len(nodes_mot)):
            if i==y:
                continue
            arr_seq=(nodes_seq[i],nodes_seq[y])
            arr_mot=(nodes_mot[i],nodes_mot[y])
            if arr_seq in sequence.edges:
                if arr_mot not in motif.edges:
                    allowed_plus=allowed_plus-len(sequence.edges[arr_seq]["label"].split('+'))
                    matched.append((arr_seq,arr_mot))
                    if allowed_plus<0:
                        return None
                else:
                    diff=len(sequence.edges[arr_seq]["label"].split('+'))-len(motif.edges[arr_mot]["label"].split('+'))
                    if diff>0:
                        allowed_plus=allowed_plus-diff
                    elif diff<0:
                        allowed_minus=allowed_minus+diff
                    if allowed_plus<0 or allowed_minus<0:
                        return None
            else:
                if arr_mot in motif.edges:
                    allowed_minus=allowed_minus-len(motif.edges[arr_mot]["label"].split('+'))
                    if allowed_minus<0:
                        return None
                    
    return matched

#fonction principale
def compare2(g,Motif,nb_miss,arr_plus,arr_moins,range_diff):
    starter=get_started(g,Motif,nb_miss,range_diff)#on creer la liste de match de taille 1
    res=[]
    while(len(starter)!=0):#tant que tout les matchs n'ont pas été traité on continue
        matched,just_added,nb,nb_range=starter.pop()#on recupere notre match
        if nb<0 or nb_range<0:#si le match a trop de penalite on l'abandonne
            continue

        if len(just_added)==0:#si tout les noeuds possible on été visite alors c'est potentiellement une solution
            if is_fully_cycled(matched):#on verifie qu'il est bien cyclique
                    
                    new_res=complete2(matched,arr_plus,arr_moins,g,Motif)#on regarde aussi que le nombre d'arete en plus et en moins respecte nos penalite
                    if new_res is not None:#si c 'est le cas on l aoute a la liste des solutions
                        res.append(new_res)
            continue
        
        else:
            graph,motif=just_added[0]#on traite la premiere arete
            a1,b1=graph#unpacking des noeuds de l'arete de la structure
            a2,b2=motif#unpacking des noeuds de l'arete du motif
            add_on=[]
            a=[(b1,i) for i in g.neighbors(b1) if check_redundancy((b1,i),matched,just_added) is False]
            c=[(b2,i) for i in Motif.neighbors(b2) if check_redundancy_motif((b2,i),matched,just_added) is False]
            permut = itertools.permutations(a, min(len(c),len(a)))#on fait la liste de tous les associations possibles des successeurs de b1 et b2
            for comb in permut:
                zipped = zip(comb, c)
                add_on.append(list(zipped))
            second_add=[]
            a=[(i,b1) for i in g.predecessors(b1) if check_redundancy((i,b1),matched,just_added) is False]
            c=[(i,b2) for i in Motif.predecessors(b2) if check_redundancy_motif((i,b2),matched,just_added) is False]
            permut = itertools.permutations(a, min(len(c),len(a)))#on fait la liste de tous les associations possibles des predecesseurs de b1 et b2
            for comb in permut:
                zipped = zip(comb, c)
                second_add.append(list(zipped))
            third_add=[]
            a=[(a1,i) for i in g.neighbors(a1) if check_redundancy((a1,i),matched,just_added) is False]
            c=[(a2,i) for i in Motif.neighbors(a2) if check_redundancy_motif((a2,i),matched,just_added) is False]
            permut = itertools.permutations(a, min(len(c),len(a)))#on fait la liste de tous les associations possibles des successeurs de a1 et a2
            for comb in permut:
                zipped = zip(comb, c)
                third_add.append(list(zipped))
            fourth_add=[]
            a=[(i,a1) for i in g.predecessors(a1) if check_redundancy((i,a1),matched,just_added) is False]
            c=[(i,a2) for i in Motif.predecessors(a2) if check_redundancy_motif((i,a2),matched,just_added) is False]  
            permut = itertools.permutations(a, min(len(c),len(a)))#on fait la liste de tous les associations possibles des successeurs de a1 et a2
            for comb in permut:
                zipped = zip(comb, c)
                fourth_add.append(list(zipped))
            
            final_add=product(product(add_on,second_add),product(third_add,fourth_add))#œn realise le produit de toutes les associations
            #final_add est donc une liste d'ensemble des agrandissmeent possible de notre match
            for i in final_add:#on parcours la liste de possibilité d'association
                nb_tmp=nb
                nb_tmp_range_diff=nb_range
                for y in i:#on va regarder pour chaque association la penalite
                    
                    a,b=y
                    cmp_missmatch,cmp_range= compare_labels(g.edges[a],Motif.edges[b])   
                    
                    if range_diff!=-1:
                        nb_tmp_range_diff=nb_tmp_range_diff-cmp_range
                    if nb_miss!=-1:
                        nb_tmp=nb_tmp-cmp_missmatch 
                if nb_tmp>=0 and nb_range>=0:#si la possibilite est toujours vraisemblable ie. elle n'a pas trop de penalité
                    copy_just_added=just_added.copy()#on copie la liste des arretes nons visité 
                    copy_just_added.remove((graph,motif))#on enleve l arete analysé a cette generation
                    copy_just_added.extend(i)#on rajoute nos possibilités trouvés à la liste d'arete non visité
                    m=matched.copy()
                    m.append((graph,motif))#l arete visite est maintenant match
                    starter.append((m,copy_just_added,nb_tmp,nb_tmp_range_diff))#on ajoute notre match à la liste
        
    return res
