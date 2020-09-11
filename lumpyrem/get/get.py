
def get_col_dict(lrout):
    with open(lrout, 'r') as file:
        a = file.readline().split()
    file.close()
    coldict = dict.fromkeys(a,None)
    return coldict

def read_cols(lrm,tuple_list):
    dictionary={}
    dictionary[lrm] = tuple_list
    return dictionary

def get_modellist(datasets):
    modellist = []
    for dset in datasets:
        modellist.append('lr_'+dset['lumprem_model_name'])