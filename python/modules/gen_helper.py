import numpy as np

QUARKS = [1, 2, 3, 4]
LEPTONS = [11, 13, 15]
NEUTRINOS = [12, 14, 16]

def status_flag(p, n):
    '''
    To save space, the status flags for each gen particle are not stored as
    individual branches, instead they are encoded into a single integer branch
    called 'GenPart_statusFlags' and can be retrieved by reading single bits
    from that number

    Example: GenPart_statusFlags == 12 (decimal) == [0...]1100 (binary)
    -> flag 0 == False
    -> flag 1 == False
    -> flag 2 == True
    -> flag 3 == True
    -> flags >=4 == False

    The list of all flags and their bit index:

    0 : isPrompt
    1 : isDecayedLeptonHadron
    2 : isTauDecayProduct
    3 : isPromptTauDecayProduct
    4 : isDirectTauDecayProduct
    5 : isDirectPromptTauDecayProduct
    6 : isDirectHadronDecayProduct
    7 : isHardProcess
    8 : fromHardProcess
    9 : isHardProcessTauDecayProduct
    10 : isDirectHardProcessTauDecayProduct
    11 : fromHardProcessBeforeFSR
    12 : isFirstCopy
    13 : isLastCopy
    14 : isLastCopyBeforeFSR

    This functions takes as input a gen particle object and the index of the status
    flag, performs a bit shift of the status flag integer to move the corresponding
    bit to the right and returns whether the resulting number is odd (last bit == 1, True)
    or even (last bit == 0, False)
    '''
    return bool((p.statusFlags >> n) % 2 != 0)

def isPrompt(p):
    return status_flag(p, 0)

def isDecayedLeptonHadron(p):
    return status_flag(p, 1)

def isTauDecayProduct(p):
    return status_flag(p, 2)

def isPromptTauDecayProduct(p):
    return status_flag(p, 3)

def isDirectTauDecayProduct(p):
    return status_flag(p, 4)

def isDirectPromptTauDecayProduct(p):
    return status_flag(p, 5)

def isDirectHadronDecayProduct(p):
    return status_flag(p, 6)

def isHardProcess(p):
    return status_flag(p, 7)

def fromHardProcess(p):
    return status_flag(p, 8)

def isHardProcessTauDecayProduct(p):
    return status_flag(p, 9)

def isDirectHardProcessTauDecayProduct(p):
    return status_flag(p, 10)

def fromHardProcessBeforeFSR(p):
    return status_flag(p, 11)

def isFirstCopy(p):
    return status_flag(p, 12)

def isLastCopy(p):
    return status_flag(p, 13)

def isLastCopyBeforeFSR(p):
    return status_flag(p, 14)

def getMassFromPDG(p):
    if abs(p.pdgId) == 12 or abs(p.pdgId) == 14 or abs(p.pdgId) == 16:
        m = 0
    elif abs(p.pdgId) == 11:
        m = 0.00051
    elif abs(p.pdgId) == 13:
        m = 0.10566
    elif abs(p.pdgId) == 15:
        m = 1.77686
    elif abs(p.pdgId) == 5:
        m = 4.18
    else:
        m = -1
    return m


def getChargeFromPDG(p):
    charge = 0.

    if abs(p.pdgId) == 2 or abs(p.pdgId) == 4 or abs(p.pdgId) == 6: # up type quarks
        charge = numpy.sign(p.pdgId) * 0.666666
    elif abs(p.pdgId) == 1 or abs(p.pdgId) == 3 or abs(p.pdgId) == 5: # down type quarks
        charge = numpy.sign(p.pdgId) * -0.333333
    elif abs(p.pdgId) == 11 or abs(p.pdgId) == 13 or abs(p.pdgId) == 15: # charged leptons
        charge = -numpy.sign(p.pdgId)
    elif abs(p.pdgId) == 24: # W boson
        charge = numpy.sign(p.pdgId)

    return charge


# checks whether in the production chain of a particle (p) up to the incoming particles a particle with pdgId in ids exists
def motherchainContains(gen_particles, p, ids):
    mother = gen_particles[p.genPartIdxMother]
    while(mother.status != 21):
        if (mother.pdgId in ids):
            return True
        mother = gen_particles[mother.genPartIdxMother]
    return False


def findGeneratorTopIdx(p, gen_particles):
    ''' find first top quark in a particle's mother decay chain '''
    motherIdx = p.genPartIdxMother
    while (motherIdx>=0):
        if abs(gen_particles[motherIdx].pdgId)==6:
            return motherIdx
        motherIdx = gen_particles[motherIdx].genPartIdxMother
    return -1
    
    
def checkMother(p, target_pdg, gen_particles):
    ''' check mother particle above p is abs(pdgId) == abs(pdg) or not '''
    mother_idx = p.genPartIdxMother
    mother_pdg = gen_particles[mother_idx].pdgId
    particle_pdg = p.pdgId

    # goes back up decay chain until it is another genparticle
    while mother_pdg == particle_pdg:
        mother_idx = gen_particles[mother_idx].genPartIdxMother
        mother_pdg = gen_particles[mother_idx].pdgId

    # check if the other genparticle has the target pdg (e.g., 24 for W boson)
    # kwargs['same_sign']
    if mother_pdg == target_pdg:
        return True, mother_pdg
    else:
        return False, mother_pdg
    
def findLastCopy(p, pdg, depth, genParticles):
    ''' recursive function to find the last copy '''
    isLastCopy_flag = False
    
    if isLastCopy(p):
        isLastCopy_flag = True
        return p, depth, isLastCopy_flag
        
    remaining_genParticles = map(lambda part_i: part_i, genParticles)[p._index+1:]
    remaining_pdg_list_np = np.array([part_i.pdgId for part_i in remaining_genParticles])
    remaining_mother_list_np = np.array([part_i.genPartIdxMother for part_i in remaining_genParticles])
    remaining_index_list_np = np.array([part_i._index for part_i in remaining_genParticles])
    
    # 1. check if any daughters
    # 2. if not then last copy
    # 3. if there are, then check pdg id and hard process
    # 4. if 0, last copy
    # 5. if 1, repeat function for that one with that p._index
    # 6. if > 1, repeat function for all of the potential daughters and choose one with maximum depth or if last copy flag
    any_daughters = np.where(remaining_mother_list_np == p._index)[0]
    
    if len(any_daughters) == 0:
        return p, depth, isLastCopy_flag
    
    # same pdg and fromHardProcess requirement
    any_daughters_pdg = remaining_pdg_list_np[any_daughters] == pdg
    any_daughters_list = remaining_index_list_np[any_daughters][any_daughters_pdg]
    any_daughters_list_hardprocess = [ fromHardProcess( genParticles[p_idx] ) for p_idx in any_daughters_list ]
    
    # if no requirements met, it is the last copy
    if sum(any_daughters_list_hardprocess) == 0:
        return p, depth, isLastCopy_flag
    
    # if one meets requirement --> recursively check if this is the last copy
    elif sum(any_daughters_list_hardprocess) == 1:
        any_daughters_idx = int( any_daughters_list[any_daughters_list_hardprocess] )
        return findLastCopy( genParticles[any_daughters_idx], pdg, depth + 1, genParticles)
        
    # if > 1 meets requirement --> recursively check all possibilities and choose deepest one or the 'isLastCopy' one (priority)
    else:
        any_daughters_idx_list = any_daughters_list[any_daughters_list_hardprocess]
        
        daughter_finalCopy_idx = []
        daughter_finalCopy_depth = []
        daughter_finalCopy_isLastCopy = []
        
        # branching to handle multiple hard process copies
        for i, weird_daughter in enumerate(any_daughters_idx_list):
            finalCopy_idx, finalCopy_depth, finalCopy_isLastCopy = findLastCopy(genParticles[weird_daughter], pdg, depth + 1, genParticles)
            daughter_finalCopy_idx.append(finalCopy_idx._index)
            daughter_finalCopy_depth.append(finalCopy_depth)
            daughter_finalCopy_isLastCopy.append(finalCopy_isLastCopy)
            
        # check if any branched daughters are isLastCopy flagged
        daughter_isLastCopy = np.where( np.array(daughter_finalCopy_isLastCopy) == True )
        if len(daughter_isLastCopy[0]) == 1:
            daughter_finalCopy = daughter_finalCopy_idx[int(daughter_isLastCopy[0])]
            daughter_depth = daughter_finalCopy_depth[int(daughter_isLastCopy[0])]
            daughter_isLastCopy_flag = True
            return genParticles[daughter_finalCopy], daughter_depth, daughter_isLastCopy_flag

        # shouldn't really be possible but would choose first one
        elif len(daughter_isLastCopy[0]) > 1:
            print("multiple last copy flagged daughters")
            daughter_finalCopy = daughter_finalCopy_idx[int(daughter_isLastCopy[0][0])]
            daughter_depth = daughter_finalCopy_depth[int(daughter_isLastCopy[0][0])]
            daughter_isLastCopy_flag = True
            return genParticles[daughter_finalCopy], daughter_depth, daughter_isLastCopy_flag
            
        # if none isLastCopy flagged --> move onto 'deepest' one in the chain
        else:
            daughter_max = np.where( np.array(daughter_finalCopy_depth) == np.array(daughter_finalCopy_depth).max() )

            # if one of these, cool
            if len(daughter_max[0]) == 1:
                daughter_finalCopy = daughter_finalCopy_idx[int( daughter_max[0] )]
                daughter_depth = daughter_finalCopy_depth[int( daughter_max[0] )]
                return genParticles[daughter_finalCopy], daughter_depth, daughter_isLastCopy_flag

            # if more than one, ffs, choose first one
            else:
                daughter_finalCopy = daughter_finalCopy_idx[int( daughter_max[0][0] )]
                daughter_depth = daughter_finalCopy_depth[int( daughter_max[0][0] )]
                return genParticles[daughter_finalCopy], daughter_depth, daughter_isLastCopy_flag
                
            # can't really have len(0) here ...
            
    # end of findLastCopy() function ---------------

def find_last_copy(p, pdg, depth, gen_particles):
    '''Recursive function to find the last copy'''
    
    if is_last_copy(p):
        return p, depth, True
        
    remaining_gen_particles = map(lambda part_i: part_i, gen_particles)[p._index+1:]
    remaining_pdg_ids = np.array([part_i.pdgId for part_i in remaining_gen_particles])
    remaining_mother_idxs = np.array([part_i.genPartIdxMother for part_i in remaining_gen_particles])
    remaining_idxs = np.array([part_i._index for part_i in remaining_gen_particles])
    
    # Find indices of any daughters
    daughter_idxs = np.where(remaining_mother_idxs == p._index)[0]
    
    # If no daughters, this is the last copy
    if len(daughter_idxs) == 0:
        return p, depth, False
    
    # Check for daughters with same pdg and fromHardProcess requirement
    daughter_pdgs = remaining_pdg_ids[daughter_idxs] == pdg
    daughter_list = remaining_idxs[daughter_idxs][daughter_pdgs]
    daughter_hardprocess_list = [ from_hard_process( gen_particles[p_idx] ) for p_idx in daughter_list ]
    
    # If no hard process daughters, this is the last copy
    if sum(daughter_hardprocess_list) == 0:
        return p, depth, False
    
    # If one hard process daughter, recursively check that daughter
    elif sum(daughter_hardprocess_list) == 1:
        daughter_idx = int( daughter_list[daughter_hardprocess_list] )
        return find_last_copy( gen_particles[daughter_idx], pdg, depth + 1, gen_particles)
        
    # If multiple hard process daughters, recursively check all and choose deepest one or last copy
    else:
        hardprocess_idxs = daughter_list[daughter_hardprocess_list]
        results = [ find_last_copy(gen_particles[idx], pdg, depth + 1, gen_particles) for idx in hardprocess_idxs ]
        final_idxs, final_depths, final_is_last_copy_flags = zip(*results)
        
        # Check for last copy daughters
        last_copy_idxs = np.where( np.array(final_is_last_copy_flags) == True )
        if len(last_copy_idxs[0]) == 1:
            return gen_particles[final_idxs[last_copy_idxs[0][0]]], final_depths[last_copy_idxs[0][0]], True

        elif len(last_copy_idxs[0]) > 1:
            print("Multiple last copy flagged daughters")
            return gen_particles[final_idxs[last_copy_idxs[0][0]]], final_depths[last_copy_idxs[0][0]], True
            
        # If none last copy, return deepest one
        else:
            deepest_idxs = np.where( np.array(final_depths) == np.array(final_depths).max() )
            return gen_particles[final_idxs[deepest_idxs[0][0]]], final_depths[deepest_idxs[0][0]], False
        
def gentop_substructures_check(gentop, top_daughters, **keys):
    # return a flag given the substructures of the gentop 
    substructure_flag = ''

    # in case of top in events:
    # check what kind of decays does the top have
    # additional key: whether the top is inside a jet or not
    
    if len(top_daughters) == 3:
        if gentop.has_hadronically_decay:
            substructure_flag = 'has_hadronicTop_' + keys['flag_is_top_inside']
        elif not gentop.has_hadronically_decay:
            substructure_flag = 'has_b_plus_lepton_fromTop_' + keys['flag_is_top_inside']
        else: 
            substructure_flag = 'has_other_' + keys['flag_is_top_inside']

    elif len(top_daughters) == 2:
        daughters_pdg = list(map(lambda top_daughter: abs(top_daughter.pdgId), top_daughters))

        if gentop.has_hadronically_decay and not any(daughter == 5 for daughter in daughters_pdg):
            substructure_flag = 'has_hadronicW_fromTop_' + keys['flag_is_top_inside']
        elif not gentop.has_hadronically_decay and not any(daughter == 5 for daughter in daughters_pdg): 
            substructure_flag = 'has_leptonicW_fromTop_' + keys['flag_is_top_inside']
        elif any(daughter == 5 for daughter in daughters_pdg) and any(daughter in QUARKS for daughter in daughters_pdg): 
            substructure_flag = 'has_b_plus_quark_fromTop_' + keys['flag_is_top_inside']
        elif any(daughter == 5 for daughter in daughters_pdg) and any(daughter in LEPTONS for daughter in daughters_pdg):  
            substructure_flag = 'has_b_plus_lepton_fromTop_' + keys['flag_is_top_inside']
        elif any(daughter == 5 for daughter in daughters_pdg) and any(daughter in NEUTRINOS for daughter in daughters_pdg): 
            substructure_flag = 'has_b_fromTop_' + keys['flag_is_top_inside']
        else: 
            substructure_flag = 'has_other_' + keys['flag_is_top_inside']

    elif len(top_daughters) == 1:
        if abs(top_daughters[0].pdgId) == 5:
            substructure_flag = 'has_b_fromTop_' + keys['flag_is_top_inside']
        elif abs(top_daughters[0].pdgId) in QUARKS:
            substructure_flag = 'has_quark_fromW_fromTop_' + keys['flag_is_top_inside']
        elif abs(top_daughters[0].pdgId) in LEPTONS:
            substructure_flag = 'has_leptonicW_fromTop_' + keys['flag_is_top_inside']
        else: 
            substructure_flag = 'has_other_' + keys['flag_is_top_inside']

    else: 
        substructure_flag = 'has_noTopDaughters_' + keys['flag_is_top_inside']

    return substructure_flag

def genW_substructures_check(W_daughters, **keys):
    # return a flag given the substructures of the genW
    # N.B. the W considered in this case is not coming from top decay!
    substructure_flag = ''

    if len(W_daughters)==2:
        daughters_pdg = list(map(lambda top_daughter: abs(top_daughter.pdgId), W_daughters))
        are_W_leptonic = all(daughter in LEPTONS for daughter in daughters_pdg)
        are_W_hadronic = all(daughter in QUARKS for daughter in daughters_pdg)

        if are_W_hadronic: 
            substructure_flag = 'has_hadronicW_not_fromTop'
        elif are_W_leptonic: 
            substructure_flag = 'has_leptonicW_not_fromTop'
        else: 
            substructure_flag = 'has_other'
    
    elif len(W_daughters)==1:
        if abs(W_daughters[0].pdgId) in QUARKS:
            substructure_flag = 'has_quark_fromW_not_fromTop'
        elif abs(W_daughters[0].pdgId) in LEPTONS:
            substructure_flag = 'has_leptonicW_not_fromTop'
        else: 
            substructure_flag = 'has_other'
    else:
        substructure_flag = 'has_other'

    return substructure_flag