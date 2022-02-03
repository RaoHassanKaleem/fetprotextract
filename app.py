from io import StringIO
from Bio import SeqIO
import pandas as pd
import streamlit as st
from PIL import Image
import extractFeatures as fe
import numpy as np
import base64
from stqdm import stqdm

icon = Image.open('fav.png')
st.set_page_config(page_title='Feature Extraction App - Proteins', page_icon = icon)

def seqValidator(seq):
    allowed_chars = set('ACDEFGHIKLMNPQRSTVWXYacdefghiklmnpqrstvwxy')
    if set(seq).issubset(allowed_chars):
        return True
    return False

seq = ""
st.header("""Feature Extraction App - Proteins""")

file_ = open("./WebPic.gif", "rb")
contents = file_.read()
data_url = base64.b64encode(contents).decode("utf-8")
file_.close()
st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="Header Image">', unsafe_allow_html=True)

st.sidebar.subheader(("Input Sequence(s) (FASTA FORMAT ONLY)"))
fasta_string  = st.sidebar.text_area("Sequence Input", height=200)
            
if st.sidebar.button("SUBMIT"):
    allFVs = []
    seqs = []
    excluded = []
    if(fasta_string==""):
        st.info("Please input the sequence first.")
    else:
        st.info("Loading Data...")
    fasta_io = StringIO(fasta_string) 
    records = SeqIO.parse(fasta_io, "fasta")
    for rec in records:
        seq_id = str(rec.id)
        seq=str(rec.seq)
        if(seqValidator(seq)):
            seqs.append(seq)
        else:
            excluded.append(seq_id)
    st.info("Computing Features...")
    i = 0
    for _ in stqdm(range(len(seqs))):
        allFVs.append(fe.calcFV(seqs[i].upper()))
        i+=1
    
    if len(allFVs)!=0:
        myDf = pd.DataFrame(allFVs)
        myCSV = myDf.to_csv(index=False, header=None).encode('utf-8')
        st.download_button(label="Download Data as CSV", data=myCSV, file_name='FVs.csv', mime='text/csv')
    if len(excluded)!=0:
        st.info("Follwoing are the IDs of excluded sequences containing invalid characters")
        st.write(excluded)
    fasta_io.close()
