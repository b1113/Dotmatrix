import streamlit as st
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

def dotmatrix(f1,f2,win):
    record1=next(SeqIO.parse(f1,"fasta"))
    record2=next(SeqIO.parse(f2,"fasta"))
    width=500 # 実際に描く幅
    height=500 # 実 際に描く高さ

    seq1=record1.seq
    seq2=record2.seq

    len1=len(seq1)-win+1 # 配列1の長さ
    len2=len(seq2)-win+1 # 配列2の長さ

    image=np.zeros((height,width)) # 実 際 に 描 く 幅 ・高さの行列

    hash={}

    for x in range(len1):
            subseq1=seq1[x:x+win]
            if subseq1 not in hash:
                hash[subseq1]=[]
            hash[subseq1].append(x)


    for y in range(len2):
        sub2=seq2[y:y+win]
        py=int(y/len2*height) # y を 画 像 位置 py に
        if sub2 in hash:
            for x in hash[sub2]:
                px=int(x/len1*width) # x を 画 像 位置 px に
                image[py,px]=1 # [py,px]に ド ッ ト


    plt.imshow(image,extent=(1,len1,len2,1),cmap="Grays")
    st.pyplot(plt)

st.title("Dot matrix") # タ イ トル
file1=st.sidebar.file_uploader("Sequence file 1:")
file2=st.sidebar.file_uploader("Sequence file 2:")
win=st.sidebar.slider("Window size:",4,100,10) 

from io import StringIO
if file1 and file2:
    with StringIO(file1.getvalue().decode("utf-8")) as f1,\
         StringIO(file2.getvalue().decode("utf-8")) as f2:
        dotmatrix(f1,f2,win)