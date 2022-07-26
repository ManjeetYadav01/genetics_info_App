import streamlit as st
import matplotlib.pyplot as plt 
import matplotlib
matplotlib.use("Agg")
from Bio.Seq import Seq 
from Bio import SeqIO
from collections import Counter
import neatbio.sequtils as utils
import numpy as np 
from PIL import Image
from io import StringIO 
import pandas as pd
from PIL import Image


def delta(x,y):
    return 0 if x == y else 1

def M(seq1,seq2,i,j,k):
    return sum(delta(x,y) for x,y in zip(seq1[i:i+k],seq2[j:j+k]))


def makeMatrix(seq1,seq2,k):
    n = len(seq1)
    m = len(seq2)
    return [[M(seq1,seq2,i,j,k) for j in range(m-k+1)] for i in range(n-k+1)]


def dotplotx(seq1,seq2):
    plt.imshow(np.array(makeMatrix(seq1,seq2,1)))
    # on x-axis list all sequences of seq 2
    xt=plt.xticks(np.arange(len(list(seq2))),list(seq2))
    # on y-axis list all sequences of seq 1
    yt=plt.yticks(np.arange(len(list(seq1))),list(seq1))
    plt.show()



activity = ['Intro','DNA Sequence','DotPlot','protien sequence',"About"]
choice = st.sidebar.selectbox("Select Activity",activity)
if choice == 'Intro':
    image = Image.open('home1.png')

    st.image(image, width=400)
    st.title("GEneticsInfo App")
    
    st.subheader("Intro")
elif choice == "DNA Sequence":
    image = Image.open('dna3.png')

    st.image(image, width=600)
    st.title("GEneticsInfo App")
    st.subheader("DNA Sequence Analysis")

    seq_file = st.file_uploader("Upload FASTA File",type=["fasta","fa"])
    if seq_file is not None:
        stringio = StringIO(seq_file.getvalue().decode("utf-8"))
        for record in SeqIO.parse(stringio, 'fasta'):
            
            dna_seq=list(str(record.seq))
            
            
            details = st.radio("Details",("Description","Sequence"))
            if details == "Description":
                st.write(record.description)
            elif details == "Sequence":
                st.write(list(str(record.seq)))
                
                
            st.subheader("Nucleotide Frequency")
            dna_freq = Counter(dna_seq)
            st.write(dna_freq)
            adenine_color = st.color_picker("Adenine Color")
            thymine_color = st.color_picker("thymine Color")
            guanine_color = st.color_picker("Guanine Color")
            cytosil_color = st.color_picker("cytosil Color")
            
            if st.button("Plot Freq"):
                barlist = plt.bar(dna_freq.keys(),dna_freq.values())
                barlist[2].set_color(adenine_color)
                barlist[3].set_color(thymine_color)
                barlist[1].set_color(guanine_color)
                barlist[0].set_color(cytosil_color)
                st.set_option('deprecation.showPyplotGlobalUse', False)

                st.pyplot()
                
            st.subheader("DNA Composition")
            gc_score = utils.gc_content(str(dna_seq))
            at_score = utils.at_content(str(dna_seq))
            st.json({"GC Content":gc_score,"AT Content":at_score})

			# Nucleotide Count
            nt_count = st.text_input("Enter Nucleotide Here","Type Nucleotide Alphabet")
            st.write("Number of {} Nucleotide is ::{}".format((nt_count),str(dna_seq).count(nt_count)))

elif choice == "DotPlot":
    image = Image.open('graph.png')

    st.image(image, width=700)
    st.title("$GEneticsInfo App$")
    st.subheader("Generate Dot Plot For Two Sequences")
    seq_file1 = st.file_uploader("Upload 1st FASTA File",type=["fasta","fa"])
    seq_file2 = st.file_uploader("Upload 2nd FASTA File",type=["fasta","fa"])

    if seq_file1 and seq_file2 is not None:
                    stringio1 = StringIO(seq_file1.getvalue().decode("utf-8"))
                    for record1 in SeqIO.parse(stringio1, 'fasta'):
                        dna_seq1=list(str(record1.seq))
                        
                    stringio2 = StringIO(seq_file2.getvalue().decode("utf-8"))
                    for record2 in SeqIO.parse(stringio2, 'fasta'):
                        dna_seq2=list(str(record2.seq))
                        
			

                    details = st.radio("Details",("Description","Sequence"))
                    if details == "Description":
                        st.write(record1.description)
                        st.write("=====================")
                        st.write(record2.description)
                    elif details == "Sequence":
                        st.write(list(str(record1.seq)))
                        st.write("=====================")
                        st.write(list(str(record2.seq)))


                    cus_limit = st.number_input("Select Max number of Nucleotide",10,200,50)
                    if st.button("Dot Plot"):
                        st.write("Comparing the first {} Nucleotide of the Two Sequences".format(cus_limit))
                        dotplotx(dna_seq1[0:cus_limit],dna_seq2[0:cus_limit])
                        st.set_option('deprecation.showPyplotGlobalUse', False)

                        st.pyplot()
            
            
elif choice == "protien sequence":
    image = Image.open('protien.png')

    st.image(image, width=700)
    st.title("GEneticsInfo App")
    df = pd.DataFrame()
    st.subheader("checking two protien sequences are equal or not ")
    seq_file1 = st.file_uploader("Upload 1st FASTA File",type=["fasta","fa"])
    seq_file2 = st.file_uploader("Upload 2nd FASTA File",type=["fasta","fa"])

    if seq_file1 and seq_file2 is not None:
        
        stringio1 = StringIO(seq_file1.getvalue().decode("utf-8"))
        for record1 in SeqIO.parse(stringio1, 'fasta'):
            dna_seq1=list(str(record1.seq))
            st.subheader("Nucleotide Frequency1")
            dna_freq1 = Counter(dna_seq1)
            a=st.write(dna_freq1)
                        
        stringio2 = StringIO(seq_file2.getvalue().decode("utf-8"))
        for record2 in SeqIO.parse(stringio2, 'fasta'):
            dna_seq2=list(str(record2.seq))
            st.subheader("Nucleotide Frequency2")
            dna_freq2 = Counter(dna_seq2)
            b=st.write(dna_freq2)
        if st.button('check for match'):    
            if dna_seq1==dna_seq2:
                st.write('its a match')
                st.balloons()
            else:
                st.write('not a match')
    
    
elif choice=='About':
    st.title('GEneticsInfo App')
    st.subheader('thanks for visiting')