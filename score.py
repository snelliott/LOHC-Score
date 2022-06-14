""" Scoring System Script for LOHC dehydrogenation energies
"""

import argparse

import automol
import ptab
from autofile import io_ as io
from autorun import execute_function_in_parallel


def get_graphs(input_file, nprocs=30):
    """get dictionary of ichs with corresponding gras
       and tanimoto scores from ichs and scores
       from Logans csv file
    :type input_file: str
    :param input_file: name of csv input file
    :type nprocs: int
    :param nprocs: number of processors

    :rtype: ({str: mol_gra}, {str: float}):
        strings are inchi identifier, mol_gra is a dictionary that
        if the molecular graph in automol format, float in the tanimoto
        simarity score read out of the input file
    """
    def _csv_line_to_graph(line):
        try:
            if '"' in line.split('.')[0]:
                ich = line.split('"')[1]
                gra = automol.inchi.graph(
                    line.split('"')[1], stereo=False)
                tan = float(line.split('"')[2].split(',')[1])
            else:
                ich = line.split(',')[1]
                gra = automol.inchi.graph(
                    line.split(',')[1], stereo=False)
                tan = float(line.split(',')[2])
            return (ich, gra, tan,)
        except:
            print('bad line', line)
            return (None, None, None)

    def _csv_lines_to_graphs(lines, output_queue):
        line_info = [_csv_line_to_graph(line) for line in lines]
        output_queue.put(tuple(line_info))

    gra_dct = {}
    tan_dct = {}
    args = ()
    lines = io.read_file(input_file).splitlines()[1:]
    ich_info = execute_function_in_parallel(
        _csv_lines_to_graphs, lines, args, nprocs=nprocs)
    for (ich, gra, tan) in ich_info:
        if ich is None:
            continue
        gra_dct[ich] = gra
        tan_dct[ich] = tan
    return gra_dct, tan_dct


def descriptor_counts(gra_dct, tan_dct, nprocs=30):
    """count the number of descriptors for each InChI

    :type gra_dct: dct{string: dct}
    :param gra_dct: dictionary of keys-inchi,
        values-molecular graph in automol format
    :type tan_dct: dct{string: float}
    :param tan_dct: dictionary of keys-inchi,
        values-tanimoto similarity score
    :type nprocs: int
    :param nprocs: numbr of processors

    :rtype dct: dictionary of key-inchi
        values-tuple of descriptor counts
    """
    def _get_rank_info(gra_dct, tan_dct, ichs, output_queue):
        rank_info_lst = ()
        for ich in ichs:
            gra = gra_dct[ich]
            rank_info = (ich,)
            rank_info += (tan_dct[ich],)
            rank_info += (percent_weight_h2(gra),)
            rank_info += carbon_hybridizations(gra)
            ring_info = ring_descriptors(gra)
            bi_atms, tri_atms = ring_info[-2:]
            rank_info += ring_info[:-3]
            rank_info += one_three_interaction(
                gra, bi_atms, tri_atms)
            rank_info_lst += (rank_info,)
        output_queue.put(rank_info_lst)
    rank_dct = {
        'keys': [
            '   score', 'tanimoto', '% wt h2',
            'sp3_carb', 'sp2_carb', 'sp_carb',
            'M-no5mem', '  M-5mem', 'M-nitcont', ' M-SOcont',
            'B-no5mem', '  B-5mem', 'B-nitcont', ' B-SOcont',
            'P-no5mem', '  P-5mem', 'P-nitcont', ' P-SOcont',
            'nonCfuse', 'SObadpos', 'SO_adjrng',
            '3mem_bi+', 'one-posi', 'onethree']}

    ichs = list(gra_dct.keys())
    args = (gra_dct, tan_dct,)
    rank_info_lst = execute_function_in_parallel(
        _get_rank_info, ichs, args, nprocs=nprocs)
    for info in rank_info_lst:
        rank_dct[info[0]] = info[1:]
    return rank_dct


def percent_weight_h2(gra):
    """ calculate %wt H2
    :type gra: dct
    :param gra: molecular graph in automol format

    :rtype: float: %wt H2
    """
    def _mass(symbs):
        return sum(list(map(ptab.to_mass, symbs)))

    atm_symbs = automol.graph.atom_symbols(gra)
    de_sat_mass = _mass(atm_symbs)
    unsat = automol.graph.unsaturated_atom_keys(gra)
    h_symbs = ['H']*len(unsat)
    h2_mass = _mass(h_symbs)
    sat_mass = de_sat_mass + h2_mass
    return (1-de_sat_mass/sat_mass)*100


def carbon_hybridizations(gra):
    """ retrieve the type of atoms, for carbons break into sp3, sp2, sp based on
        unsaturation level
    :type gra: dct
    :param gra: molecular graph in automol format

    :rtype: tuple(int): tuple of number of each hybriziation
    """
    sp3_carb_atms = 0
    sp2_carb_atms = 0
    sp_carb_atms = 0
    atm_vals = automol.graph.atom_unsaturated_valences(gra)
    # unsat_atms = automol.graph.unsaturated_atom_keys(gra)
    for key, atms in automol.graph.atom_symbol_keys(
            automol.graph.explicit(gra)).items():
        for atm in atms:
            if key == 'C':
                if atm_vals[atm] == 0:
                    sp3_carb_atms += 1
                elif atm_vals[atm] == 1:
                    sp2_carb_atms += 1
                elif atm_vals[atm] == 2:
                    sp_carb_atms += 1
    return (
        sp3_carb_atms, sp2_carb_atms, sp_carb_atms)


def one_three_interaction(gra, bi_atms, tri_atms):
    """returns if there is a 1,3 interaction
       or an N,3 interaction
    """
    rings_atms = list(automol.graph.rings_atom_keys(gra))
    atm_dct = automol.graph.atoms(gra)
    ngb_dct = automol.graph.atoms_neighbor_atom_keys(gra)
    rank_one = 0
    rank_one_three = 0
    joining_bi_atms = [
        atm for n, atm in enumerate(bi_atms) if atm in bi_atms[:n]]
    joining_tri_atms = [
        atm for n, atm in enumerate(tri_atms) if atm in tri_atms[:n]]
    joint_atms = joining_bi_atms + joining_tri_atms
    for ring_atms in rings_atms:
        one = False
        for idx, atm_1 in enumerate(ring_atms[:-2]):
            atm_3 = ring_atms[idx+2]
            if any(joint_atm in [atm_1, atm_3] for joint_atm in joint_atms):
                continue
            one_three = True
            if not (atm_dct[atm_1][0] == 'N' or any(
                    (atm_dct[ngb][0] != 'H' and ngb not in ring_atms)
                    for ngb in ngb_dct[atm_1])):
                one_three = False
            elif any((
                    atm_dct[ngb][0] != 'H' and ngb not in ring_atms)
                    for ngb in ngb_dct[atm_1]):
                one = True
            if not (atm_dct[atm_3][0] == 'N' or any(
                    (atm_dct[ngb][0] != 'H' and ngb not in ring_atms)
                    for ngb in ngb_dct[atm_3])):
                one_three = False
            elif any((
                    atm_dct[ngb][0] != 'H' and ngb not in ring_atms)
                    for ngb in ngb_dct[atm_3]):
                one = True
            if one_three:
                rank_one_three += 1
        if one:
            rank_one += 1
    return rank_one, rank_one_three


def intersection(lst1, lst2):
    """returns intersec of two lists
    """
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


def _polycyclic_rings(rings_atms):
    """ returns tuple of monocyclic rings for a tuple
        of all rings in a tuple(tuple(int,)) where int is atm key
    """
    # find each ring that each ring is connected to
    connected_ring_idxs_tup = ()
    for ring_atms_i in rings_atms:
        connected_ring_idxs = set()
        for idx_j, ring_atms_j in enumerate(rings_atms):
            if any([j_atm in ring_atms_i for j_atm in ring_atms_j]):
                connected_ring_idxs.add(idx_j)
        connected_ring_idxs_tup += (connected_ring_idxs,)

    # gotta keep walking through the connections for polycyclics
    missing_conns = True
    while missing_conns:
        new_conns_tup = ()
        for conns_i in connected_ring_idxs_tup:
            new_conns = set()
            for conns_j in connected_ring_idxs_tup:
                if any([j_idx in conns_i for j_idx in conns_j]):
                    [new_conns.add(j_idx) for j_idx in conns_j]
            new_conns_tup += (new_conns,)
        if connected_ring_idxs_tup == new_conns_tup:
            missing_conns = False
        else:
            connected_ring_idxs_tup = new_conns_tup

    # remove redundant ring groups
    unique_conn_idxs = ()
    for idxs in connected_ring_idxs_tup:
        if idxs not in unique_conn_idxs:
            unique_conn_idxs += (idxs,)

    # split into mono, bi, and tri+ cyclic groups
    ret = [(), (), ()]
    for ring_network in unique_conn_idxs:
        ring_atm_network = ()
        for idx in ring_network:
            ring_atm_network += (rings_atms[idx],)
        ret[min(len(ring_network)-1, 2)] += (ring_atm_network,)
    return ret


def _sub_ring_types(ring_networks, atm_dct):
    """ break ring networks into number of:
        non 5-member ring containing network
        5-member ring containing network
        N-atm containing network
        S/O-atm containing network
    """
    num_non5 = 0
    num_5 = 0
    num_Natm = 0
    num_SOatm = 0
    for ring_network in ring_networks:
        if any([len(ring) == 5 for ring in ring_network]):
            num_5 += 1
        else:
            num_non5 += 1
        all_rings_atms = sum(ring_network, ())
        if any([atm_dct[atm][0] == 'N' for atm in all_rings_atms]):
            num_Natm += 1
        if any([atm_dct[atm][0] in ['S', 'O'] for atm in all_rings_atms]):
            num_SOatm += 1
    return (num_non5, num_5, num_Natm, num_SOatm)


def _fuse_points(ring_networks):
    """ atm keys at which rings are fused
    """
    fuse_keys = ()
    for ring_network in ring_networks:
        for idx, ring in enumerate(ring_network):
            other_ring_atms = sum(ring_network[:idx] + ring_network[idx+1:], ())
            ring_fuse_points = [atm for atm in ring if atm in other_ring_atms]
            fuse_keys += tuple(ring_fuse_points)
    return tuple(set(fuse_keys))



def ring_descriptors(gra):
    """ heterocyclic rank: number of overlapping rings
    """
    rings_atms = list(automol.graph.rings_atom_keys(gra))
    atm_dct = automol.graph.atoms(gra)
    ngb_dct = automol.graph.atoms_neighbor_atom_keys(gra)

    # split up rings into monocyclic, bicyclic, and polycyclic (which is >2)
    monocyc_rings, bicyc_rings, polycyc_rings = _polycyclic_rings(rings_atms)
    mono_atms = tuple(set(sum(sum(monocyc_rings, ()), ())))
    bi_atms = tuple(set(sum(sum(bicyc_rings, ()), ())))
    poly_atms = tuple(set(sum(sum(polycyc_rings, ()), ())))

    bifuse_keys = _fuse_points(bicyc_rings)
    polyfuse_keys = _fuse_points(polycyc_rings)
    # strain where rings fuse if its not a carbon atom
    num_noncarb_fuse_atms = 0
    for key in bifuse_keys + polyfuse_keys:
        if atm_dct[key][0] != 'C':
            num_noncarb_fuse_atms += 1

    # if S or O is in bi+cyclics, it should be adjacent to a fuse point
    num_so_badpos = 0
    for key in bi_atms + poly_atms:
        if atm_dct[key][0] in ['S', 'O']:
            if not any([
                    ngb_atm in bifuse_keys + polyfuse_keys
                    for ngb_atm in ngb_dct[key]]):
                num_so_badpos += 1

    # S or O causes ring strain if present in two connected five membered rings
    num_so_adj_rings = 0
    for fuse_key in bifuse_keys + polyfuse_keys:
        fused_rings = [
            ring_atms for ring_atms in rings_atms if fuse_key in ring_atms]
        adjso = 0
        for fused_ring in fused_rings:
            if any([atm_dct[ring_atm] in ['S', 'O'] for ring_atm in fused_ring]):
                adjso += 1
        if adjso > 1:
            num_so_adj_rings += 1

    # number of 3-member rings in bi+cyclics
    num_3mem_nonmono = 0
    threemem_rings = [ring_atms for ring_atms in rings_atms if len(
        ring_atms) == 3]
    for ring_atms in threemem_rings:
        if any([ring_atm in bi_atms + poly_atms for ring_atm in ring_atms]):
            num_3mem_nonmono += 1

    return (
        _sub_ring_types(monocyc_rings, atm_dct) +
        _sub_ring_types(bicyc_rings, atm_dct) +
        _sub_ring_types(polycyc_rings, atm_dct) +
        (num_noncarb_fuse_atms,) +
        (num_so_badpos, num_so_adj_rings,) +
        (num_3mem_nonmono,) +
        (mono_atms, bi_atms, poly_atms,))


def score(rank_dct):
    """ sort the inchis in the rank_dct dictionary according to a score
        assigned according to fitted weight factors for each group in the
        rank_dct
    :type rank_dct: dct
    :param rank_dct: dictionary with key-ich and values-list of
        counts for each descriptor

    :rtype list(list(float, )): a list of lists of descriptor
        counts sorted by the score of that descriptor count set
    """
    factors = [
        # tanimoto, % wt h2,
        0, 0,
        # sp3_carb, sp2_carb, sp_carb,
        0.1153, 0.2872, 1.1298,
        # M-no5mem, M-5mem, M-nitcont,  M-SOcont,
        -1.2127, -1.6513, -0.0955, 0.9387,
        # B-no5mem, B-5mem, B-nitcont,  B-SOcont,
        5.2805, 5.0358, -0.3393, -0.0191,
        # P-no5mem, P-5mem, P-nitcont,  P-SOcont,
        1.0000, 3.4910, 0.5179, 0.6033,
        # nonCfuse, SObadpos, SO_adjrng,
        0.5200, 0.9095, 2.0000,
        # 3mem_bi+, one-posi, onethree]}
        3.9448, -0.0735, -0.2623]
    rank_lst = []
    for ich, rank_info in rank_dct.items():
        if ich == 'keys':
            continue
        # factor is 0 if %wt h2 is too low
        if rank_info[1] > 5.4:
            score_ = sum([i * j for i, j in zip(rank_info, factors)])
        else:
            score_ = 0.
        rank_lst.append([ich, score_, *rank_info])
    sorted_rank_lst = sorted(
        rank_lst, key=lambda x: ((x[1]), x[2]), reverse=True)
    return sorted_rank_lst


def _command_line(ranked_lst, rank_dct):
    """ print to command line instead of to an output file
    """
    print('InChI\t\t\t\t,', ','.join(rank_dct['keys']))
    for info in ranked_lst:
        print(
            '"' + info[0] + '"' +
            ',{:.3f},'.format(info[1]) + ','.join(
                ['{:g}'.format(info) for info in info[2:]]))


def write_output(ranked_lst, rank_dct, out_file):
    """ write scores and features to an output file in csv format
        'InChI,score,tanimoto,%wtH2,...'

    :type ranked_lst: list(list(float, ))
    :param ranked_lst: a list of lists of descriptor
        counts sorted by the score of that descriptor count set
    :type rank_dct: dct
    :param rank_dct: dictionary with key-ich and values-list of
        counts for each descriptor
    :type output_file: str
    :param output_file: name for csv output file
    """

    if out_file == 'command_line':
        _command_line(ranked_lst, rank_dct)
    else:
        out_str = 'InChI\t\t\t\t,' + ','.join(rank_dct['keys'])
        out_str += '\n'
        for info in ranked_lst:
            out_str += (
                '"' + info[0] + '"' +
                ',{:.3f},'.format(info[1]) + ','.join(
                  ['{:g}'.format(info) for info in info[2:]]) + '\n')

        io.write_file(out_file, out_str)


def main(input_file, output_file):
    """scores each molecule according to the count of each descriptor it has
       and a predetermined weights for those descriptors.

    :type input_file: string
    :param input_file: name of a csv file that contains
        <InChI>,<tanimoto similarity>
    """
    gra_dct, tan_dct = get_graphs(input_file, nprocs=30)
    rank_dct = descriptor_counts(gra_dct, tan_dct, nprocs=30)
    ranked_lst = score(rank_dct)
    write_output(ranked_lst, rank_dct, output_file)


def get_args():
    """ argparse
    """
    parser = argparse.ArgumentParser(description='Score LOHC')
    parser.add_argument(
        '-i', '--input_file', type=str,
        help='name of csv input file',
        default='input.csv')
    parser.add_argument(
        '-o', '--output_file', type=str,
        help='name of csv output file',
        default='output.csv')
    return parser.parse_args()


if __name__ == '__main__':
    ARGS = get_args()
    main(ARGS.input_file, ARGS.output_file)
