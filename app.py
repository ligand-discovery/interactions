# streamlit_app.py
import streamlit as st
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
from io import BytesIO
import xlsxwriter
import os
import sys

# relative imports
ROOT = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(ROOT, "./src/"))
from agstyler import PINLEFT, PRECISION_TWO, draw_grid

st.set_page_config(
    page_title="Ligand-Protein interactions in Chemical Proteomics Screening", 
    page_icon=":home:",
    layout="wide", # "centered",
    initial_sidebar_state="expanded"
)

st.markdown("""
  <style>
    .css-13sdm1b.e16nr0p33 {
      margin-top: -75px;
    }
  </style>
""", unsafe_allow_html=True)

hide_streamlit_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            header {visibility: hidden;}
            </style>
            """
st.markdown(hide_streamlit_style, unsafe_allow_html=True) 

pIdDf = pd.read_csv(os.path.join(ROOT, "./data/proteinNames.tsv"), sep="\t")

pId = pIdDf['UniProtID'].values
pIdDes = pIdDf['Description'].values

def applyFilters(df, pFil, pAdjFil, hitFil):
  if pFil != 'no filter':
    if pFil == '< 0.05':
      df = df[df['ml10p'] > 1.30103]
    else:
      df = df[df['ml10p'] > 2]
      
  if pAdjFil != 'no filter':
    if pAdjFil == '< 0.05':
      df = df[df['ml10adjP'] > 1.30103]
    elif pAdjFil == '< 0.1':
      df = df[df['ml10adjP'] > 1]
    else:
      df = df[df['ml10adjP'] > 0.60206]
      
  if hitFil != 'no filter':
    if hitFil == 'Low':
      df = df[df['mdfClass'] >= 1]
    elif hitFil == 'Medium (hits)':
      df = df[df['mdfClass'] >= 2]
    elif hitFil == 'Low (hits)':
      df = df[df['mdfClass'] >= 1]
    else:
      df = df[df['mdfClass'] == 3]
      
  return df

def getVarText(df):
  if (len(df.index)) > 0:
    bestProt = df["geneName"].values[0]
    numProtHitss = len(df.index)
    df.index = np.arange(1,len(df)+1)
    protList = df.index[df["accession"]==myPid].tolist()
    if len(protList) > 0:
      protRank = protList[0]
      varText1 = "hit rank is"
      varText2 = "is best"
    else:
      varText1 = "is not a hit"
      protRank = ""
      varText2 = "protein is best"
    del protList
  else:
    bestProt = "No "
    numProtHitss = 0
    varText1 = "is not a hit"
    protRank = ""
    varText2 = ""
  return [bestProt, numProtHitss, protRank, varText1, varText2]


st.sidebar.title("Ligand-Protein Interactions")
st.title("Chemical Proteomics Screening")

help_input3='''
    
    use **:blue[UniProt Accession]**, Short Gene Name(s) or Protein Description to search\n
    **Tip**:\n
    To change selected protein, **:red[NO]** need to select whole existing term, delete and type new.\n
    :blue[Just start to type new protein, old text will be automatically cleared]'''

pIdIndex = st.sidebar.selectbox(label = "Select Protein", help = help_input3, options = range(len(pIdDes)), format_func= lambda x: pIdDes[x], index= 7109)

myPid = pId[pIdIndex]

intDfOri = pd.read_csv(os.path.join(ROOT, "./data/finalScreen.tsv"), sep="\t")
fpDf = pd.read_csv(os.path.join(ROOT, "./data/finalFp.tsv"), sep="\t")
intDf = intDfOri[intDfOri["accession"]==myPid]

if len(intDf) == 0:
  st.sidebar.write("We did **:red[not]** detect selected protein interacting with any ligand in our screen, try another protein")
else:
  selectedGeneName = intDf["geneName"].values[0]
  tempDf = applyFilters(intDf, '< 0.05', '< 0.25', 'Medium (hits)')
  if (len(tempDf.index)) > 0:
    tempDf = tempDf.sort_values(by=['protHits', 'l2fc'], ascending=[True, False])
    bestFrag = tempDf["fragId"].values[0]
    top5Frags = tempDf["fragId"].values[0:5]
    numLigaHits = len(tempDf.index) # numLigaHits is already present in base input table
    varText3 = "is best"
  else:
    bestFrag = "No"
    varText3 = ""
    numLigaHits = 0
    
######## Screening Protein Centric View ############      

  st.write("**Selected Protein**: ", pIdDes[pIdIndex])
  st.markdown("""---""")

  st.subheader(f"First generation ligands (Gen1) that enrich **:blue[{selectedGeneName}]** over background")
  
  numInt = len(intDf.index)
  st.write(f"**:blue[{numInt}]** (out of 407 screened) Gen1 ligands enrich **{selectedGeneName}**. **:blue[{numLigaHits}]**/{numInt} ligands are labelled as **hits** by applying **medium** filter Set **(:blue[fS])**. **:blue[{bestFrag}]** {varText3} **hit**.")

  if (numLigaHits/407)>0.1:
    hitRatio = np.round((numLigaHits/407)*100, 1) 
    st.write(f"**:blue[{selectedGeneName}]** is a **:red[promiscuous]** protein (**hit**/enriched ratio is **:red[{hitRatio}]**%).")

  col1, col2, col3, colX, colY, colZ = st.columns(6)
  with col1:
    pFilter = st.selectbox(label = "*P* Value", help = "Select threshold for signifiance", options = ('< 0.05', 'no filter', '< 0.01'))
  with col2:
    pAdjFilter = st.selectbox(label = "adjusted *P* Value", help = "Select threshold for signifiance", options = ('< 0.25', '< 0.1', 'no filter', '< 0.05'))
  with col3:
    help_input='''
    
    **:blue[0]**.  no filter\n
    **:blue[1]**.  Low Confidence: Fc > 1,   Median > 1,   p < 0.05, adj.p < 0.25, Rank < 500\n
    **:blue[2]**.  Medium confidence ('**:blue[hits]**'): Fc > 2.3, Median > 1,   p < 0.05, adj.p < 0.25, Rank < 500\n
    **:blue[3]**.  High Confidence (also '**:blue[hits]**'):  Fc > 2.3, Median > 2.3, p < 0.01, adj.p < 0.1,  Rank < 500'''
    mdfClass = st.selectbox(label = "filter Set (**:blue[fS]**)", help = help_input, options = ('Medium (hits)', 'no filter', 'Low', 'High (hits)'))
    
  if len(tempDf.index) == 0:
    st.write("**:red[No]** data to display with selected filters. Applied **:blue[no filter]**")
    intDf = applyFilters(intDf, 'no filter', 'no filter', 'no filter')
    
  else:
    intDf = applyFilters(intDf, pFilter, pAdjFilter, mdfClass)

  del tempDf
  
  intDf = intDf.sort_values(by=['protHits', 'l2fc'], ascending=[True, False])
    
  col4, col5  = st.columns(2)
  with col4:
    formatter = {
      'fragId': ('Ligand', {**PINLEFT, 'width': 10}),
      'l2fc': ('Fc(log2)', {**PRECISION_TWO, 'width': 15}),
      'l2fcM': ('Fc Median adjusted', {**PRECISION_TWO, 'width': 25}),
      'protHits': ('# Protein Hits', {'width': 15}),
      'mdfClass': ('fS', {'width': 10})
    }
    data = draw_grid(intDf, formatter=formatter, fit_columns=True, selection='none', max_height=340)
  with col5:
    st.image(os.path.join(ROOT, "./assets/proteinCentric/") + myPid + ".png")
    
  fragId = st.sidebar.selectbox(label = "Select Gen1 Ligand", options = intDf["fragId"])
  intDf2 = intDfOri[intDfOri["fragId"]==fragId]
  
############ Screening Ligand Centric ###############################
  
  st.subheader(f"Proteins enriched by **:blue[{fragId}]**")
  
  tempDf2 = intDfOri[intDfOri["fragId"]==fragId]
  numProtDetected = len(tempDf2.index)
  tempDf2 = applyFilters(tempDf2, '< 0.05', '< 0.25', 'Medium (hits)')

  tempDf2 = tempDf2.sort_values(by=['ligHits', 'l2fc'], ascending=[True, False])
  [bestProt, numProtHits, protRank, varText, varText2] = getVarText(tempDf2)
  
  if len(tempDf2.index) == 0:
    intDf3 = applyFilters(intDf2, 'no filter', 'no filter', 'no filter')
    
  else:
    intDf3 = applyFilters(intDf2, pFilter, pAdjFilter, mdfClass)
  
  intDf3 = intDf3.sort_values(by=['ligHits', 'l2fc'], ascending=[True, False])
    
  st.sidebar.image(os.path.join(ROOT, "../assets/fragFiguresSingle/") + fragId + ".png")
  
  st.write(f"**:blue[{numProtDetected}]** proteins were enriched by ligand **{fragId}** (Fc compared to **CRF** control). **:blue[{numProtHits}]** of those proteins were labelled as **hits** by applying **medium** filter Set **(:blue[fS])**. **:blue[{bestProt}]** {varText2} **hit**. **:blue[{selectedGeneName}]** {varText} **:blue[{protRank}]**.")
  
  if (numProtHits/numProtDetected)>0.05:
    fragHitRatio = np.round((numProtHits/numProtDetected)*100, 1) 
    st.write(f"**:blue[{fragId}]** is **:red[promiscuous]** ligand (**hit**/enriched ratio is **:red[{fragHitRatio}]**%).")
  
  col6, col7  = st.columns(2)
  with col6:
    st.image(os.path.join(ROOT, "./assets/ligandVolcanoPlots/") + fragId + ".png")
  with col7:
    if len(tempDf2.index) == 0:
      st.write("**:red[No]** data to display with selected filters. Applied **:blue[no filter]**")
    formatter = {
      'accession': ('Protein', {**PINLEFT, 'width': 15}),
      'geneName': ('Gene', {**PINLEFT, 'width': 15}),
      'l2fc': ('Fc(log2)', {**PRECISION_TWO, 'width': 15}),
      'l2fcM': ('Fc Median adjusted', {**PRECISION_TWO, 'width': 25}),
      'ligHits': ('# Ligand Hits', {'width': 15}),
      'mdfClass': ('fS', {'width': 10})
    }
    data = draw_grid(
      intDf3, formatter=formatter, fit_columns=True, selection='none', max_height=340)
  if not isinstance(protRank, str):
    if protRank < 5:
      st.subheader(f"**:blue[{fragId}-{selectedGeneName}]** interaction: :first_place_medal:")
      st.write(f"**:blue[{fragId}]** is in top 5 **Ligand hits** for **{selectedGeneName}**. **:blue[{selectedGeneName}]** is in top 5 **Protein hits** for **{fragId}**.")

  del tempDf2
############# Fingerprinting / Elaborates Data ######################

  gen2List = ["C027", "C028", "C044", "C046", "C064", "C115", "C127", "C160", "C179", "C186", "C197", "C219", "C240", "C270", "C275", "C303", "C310", "C320", "C378", "C391"]
  if fragId in gen2List:
    st.markdown("""---""")
    # st.sidebar.markdown("""---""")
    
############ Elaborates Protein Centric View ##########################
    st.subheader(f"Second generation ligands (Gen2) of **:blue[{fragId}]** that compete **:blue[{selectedGeneName}]**")

    gen1Df = fpDf[fpDf["gen1Lig"]==fragId]
    
    numGen2Ligs = len(np.unique(gen1Df['fragId']))
    
    temp4Df = gen1Df[gen1Df["accession"]==myPid]
    temp4Df = temp4Df[temp4Df["mdfClass"]>=1]
    # temp4Df = temp4Df.sort_values(by='l2fc', ascending=True)
    temp4Df = temp4Df.sort_values(by='hitRank', ascending=True)
    
    sidebarList1 = temp4Df["fragId"]

    if len(temp4Df.index)>0:
      bestGen2Lig = temp4Df['fragId'].values[0]
      varText5 = "is best Gen2 ligand hit."
    else:
      bestGen2Lig = ""
      varText5 = ""
    
    st.write(f"**:blue[{numGen2Ligs}]** Gen2  ligands were screened in **competition** experiments against Gen1 ligand **{fragId}**. **:blue[{len(temp4Df.index)}]**/{numGen2Ligs} Gen2 ligands of **{fragId}** pass **low** filter Set (**:blue[fS2]**).")
    # **:blue[{bestGen2Lig}]** {varText5}")
    # **compete** **:blue[{selectedGeneName}]** after applying 
    
    formatter = {
      'fragId': ('Gen2', {**PINLEFT, 'width': 10}),
      'l2fc': ('Fc(log2)', {**PRECISION_TWO, 'width': 10}),
      'l2fcM': ('Fc Median adjusted', {**PRECISION_TWO, 'width': 20}),
      'protHits': ('# Gen2 Protein Hits', {'width': 15}),
      'mdfClass': ('fS2', {'width': 10})
    }
    
    col10, col11 = st.columns(2)
    with col10:
      st.write(f":blue[Hits] (fS2 > 0)")
      if len(temp4Df.index)>0:
        data = draw_grid(
          temp4Df, formatter=formatter, fit_columns=True, selection='none')
      
    temp4Df = gen1Df[gen1Df["accession"]==myPid]
    temp4Df = temp4Df[temp4Df["mdfClass"] < 1]
    temp4Df = temp4Df.sort_values(by='l2fc', ascending=True)
  
    sidebarList2 = temp4Df["fragId"]
    
    with col11:
      st.write(f":orange[not] Hits (fS2 = 0)")
      data = draw_grid(
        temp4Df, formatter=formatter, fit_columns=True, selection='none')
    
    temp4Df = gen1Df[gen1Df["accession"]==myPid]
    temp4Df = temp4Df.sort_values(by='l2fc', ascending=True)
    temp4Df = temp4Df.sort_values(by=['mdfClass', 'l2fc'], ascending=[False, True])

    sideBarList = pd.concat([sidebarList1, sidebarList2], sort=False)
    
############ Elaborates Side Bar Selection ##########################      

    # gen2Id = st.sidebar.selectbox(label = "Select Gen2 Ligand", options = temp4Df["fragId"])
    gen2Id = st.sidebar.selectbox(label = "Select Gen2 Ligand", options = sideBarList)
    st.sidebar.image(os.path.join(ROOT, "./assets/fragFiguresSingle/") + gen2Id + ".png")
    
############ Elaborates Ligand Centric View ##########################      

    st.subheader(f"Proteins competed by **:blue[{gen2Id}]**")
    
    gen2Df = gen1Df[gen1Df["fragId"]==gen2Id]
    
    tempDf3 = applyFilters(gen2Df, '< 0.05', '< 0.25', 'Low')
    tempDf3 = tempDf3.sort_values(by=['ligHits', 'l2fc'], ascending=[True, True])
    
    [bestProt2, numProtHits2, protRank2, varText3, varText4] = getVarText(tempDf3)
    
    st.write(f"**:blue[{len(gen2Df.index)}]** proteins were reduced in **:blue[{gen2Id}] competition** experiment (Fc compared to **{fragId}** control). **:blue[{numProtHits2}]** of those proteins are labelled as **hits** by applying **low** filter Set **(:blue[fS2])**. **:blue[{bestProt2}]** {varText4} **hit**. **:blue[{selectedGeneName}]** {varText3} **:blue[{protRank2}]**.")
    
    col1, col2, col3, colX, colY, colZ = st.columns(6)
    with col1:
      pFilterFP = st.selectbox(label = "*P* Value", help = "Select threshold for signifiance", options = ('< 0.05', 'no filter', '< 0.01'), key = 'pFilterFP')
    with col2:
      pAdjFilterFP = st.selectbox(label = "adjusted *P* Value", help = "Select threshold for signifiance", options = ('< 0.25', 'no filter', '< 0.1', '< 0.05'), key = 'pAdjFilterFP')
    with col3:
      help_input2='''
      
      **:blue[0]**.  no filter\n
      **:blue[1]**.  Low Confidence ('**:blue[hits]**'): Fc < -1,   p < 0.05, adj.p < 0.25, Rank < 500\n
      **:blue[2]**.  Medium confidence (also '**:blue[hits]**'): Fc < -1.65, p < 0.05, adj.p < 0.25, Rank < 500\n
      **:blue[3]**.  High Confidence (also '**:blue[hits]**'):  Fc < -2.3, p < 0.01, adj.p < 0.1,  Rank < 500'''
      mdfClassFP = st.selectbox(label = "Gen2 Ligand filter Set (**:blue[fS2]**)", help = help_input2, options = ('Low (hits)', 'Medium (hits)', 'no filter', 'High (hits)'), key = 'mdfClassFP')
    
    if len(tempDf3.index) == 0:
      st.write("**:red[No]** data to display with selected filters. Applied **:blue[no filter]**")
      gen2Df2 = applyFilters(gen2Df, 'no filter', 'no filter', 'no filter')
      
    else:  
      gen2Df2 = applyFilters(gen2Df, pFilterFP, pAdjFilterFP, mdfClassFP)
      
    gen2Df2 = gen2Df2.sort_values(by=['ligHits', 'l2fc'], ascending=[True, True])
    
    col8, col9  = st.columns(2)
    with col8:
      formatter = {
      'accession': ('Protein', {**PINLEFT, 'width': 15}),
      'geneName': ('Gene', {**PINLEFT, 'width': 15}),
      'l2fc': ('Fc(log2)', {**PRECISION_TWO, 'width': 10}),
      'l2fcM': ('Fc Median adjusted', {**PRECISION_TWO, 'width': 20}),
      'ligHits': ('# Gen2 Ligand Hits', {'width': 15}),
      'mdfClass': ('fS2', {'width': 10})
      }
      data = draw_grid(
        gen2Df2, formatter=formatter, fit_columns=True, selection='none', max_height=340)
    with col9:
      st.image(os.path.join(ROOT, "./assets/gen2VolcanoPlots/") + gen2Id + ".png")