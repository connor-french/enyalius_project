# create an imap for running PCAs in ipyrad

def create_imap(splist, locs_df):
  id_df = locs_df[["species", "id_code"]]
  id_df = id_df[id_df["species"].isin(splist)]
  
  imap = {}
  
  for i in splist:
    imap[i] = id_df["id_code"][id_df["species"] == i].tolist()
    
  return imap
