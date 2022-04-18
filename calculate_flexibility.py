###########################################
######## CALCULATE FLEXIBILITY ############
###########################################

"""B-values calculated for each residue"""

def calculate_flex_pos(B_val_dict, pos_dict, score_dict):

    '''Function that takes a dictionary with normalized B-values, a dictionary 
    with the indexed position of a MSA and a dictionary with BLAST scores and returns
    a list of (position, residue, flexibility, confidence) tuples for a query.'''

    flexibility_position_list = []

    query_position = pos_dict.pop("query")

    pos = 1 # counter for position

    for query_index_tuple in query_position:
        
        B_vals_weighted = []
        scores = []
        
        query_residue = query_index_tuple[0]
        query_MSA_index = query_index_tuple[1]

        for template_ID, template_index_list in pos_dict.items():
            for template_index_tuple in template_index_list:

                template_MSA_index = template_index_tuple[1]
                template_PDB_index = template_index_tuple[2]

                if query_MSA_index == template_MSA_index:
                    B_val_raw = float(B_val_dict[template_ID]["B_val_list"][template_PDB_index][1])
                    B_vals_weighted.append(B_val_raw * float(score_dict[template_ID]))
                    scores.append(float(score_dict[template_ID]))
                       
        if (B_vals_weighted and scores):
            flexibility_position_list.append((pos, query_residue,(sum(B_vals_weighted)/sum(scores)), len(scores)))

        else:
            flexibility_position_list.append((pos, query_residue,"?", 0))
        
        pos+=1
                
    return(flexibility_position_list)

    
        
def complete_missing_scores (incomplete_flexibility_list):

    '''Function that gives a flexibility score to regions without homology. It uses the mean
    flexibility score for each residue in the query sequence.'''

    # complete list

    complete_flexibility_list = []

    # calculate mean flexibilities for each amino acid

    residues_list = ["V", "N", "Q", "L", "P", "F", "I", "Y", "W", "S", "M", "C", "T", "G", "A", "D", "E", "K", "R", "H", "X"]
    mean_flex_dict = {}

    for res in residues_list:

        scores = []

        for flex_tuple in incomplete_flexibility_list:
            if (flex_tuple[1] == res and flex_tuple[2] != "?"):
                scores.append(flex_tuple[2])

        if scores: # avoid zero division error when there is an amino acid missing
            mean_flex_dict[res]= sum(scores)/len(scores)

    for flex_tuple in incomplete_flexibility_list:
        if flex_tuple[2] == "?":
            complete_flexibility_list.append((flex_tuple[0], flex_tuple[1], mean_flex_dict.setdefault(flex_tuple[1], 0), 0))

        else:
            complete_flexibility_list.append(flex_tuple)

    return(complete_flexibility_list)
