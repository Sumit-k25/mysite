from modlamp.descriptors import PeptideDescriptor, GlobalDescriptor
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from datetime import datetime, timezone, timedelta
import os, base64, requests, mpld3, numpy as np
from modlamp.plot import helical_wheel
from Bio.Blast import NCBIXML
from datetime import datetime
from Bio.SeqUtils import seq3
from Bio.Blast import NCBIWWW
import Bio.SeqUtils.ProtParam
import modlamp.plot as pl
from Bio import SeqIO
import collections
import random
import auth
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker



def updateVisitorCount(userCount:int):
    try:
        userCount = userCount + 1
        with open("mysite/userCount.txt", "w") as file:
            file.write(str(userCount))
        return userCount
    except:
        return None



def visitorCountInit():
    try:
        userCount=0
        if os.path.exists("mysite/userCount.txt"):
            with open("mysite/userCount.txt","r") as file:
                content=file.read()
                userCount=int(content)
        else:
            with open("mysite/userCount.txt","w") as file:
                file.write("0")
                userCount=0
        return userCount
    except:
        return None



def uploadImg(file_name:str) -> str:
    link=None
    try:
        with open(file_name,'rb') as file:
            payload={
                    "key": auth.imgbb_API_Key,
                    "image": base64.b64encode(file.read()),
                    }
            response = requests.post("https://api.imgbb.com/1/upload", payload)

        # Get the URL of the uploaded image from the API response and update the 'glist' dictionary
        link = response.json()['data']['url']
        return link

    except Exception as e:
        print(f"\n\nError in uploadImg(): {e}\n\n")



def uploadTxt(file_name:str) -> str:
    link=None
    try:
        # Set up your Pinata API credentials
        api_key = auth.ipfs_API_Key
        api_secret_key = auth.ipfs_API_Secret

        # Set up the expiration time for your pin
        expiration_time = datetime.now(timezone.utc) + timedelta(hours=0.01)

        # Upload your file to Pinata and set an expiration time for the pin
        with open(file_name, 'rb') as f:
            headers = {
                        'pinata_api_key': api_key,
                        'pinata_secret_api_key': api_secret_key,
                      }
            response = requests.post(
                'https://api.pinata.cloud/pinning/pinFileToIPFS',
                headers=headers,
                files={'file': f},
                json={
                    'pinataMetadata': {
                        'expiration': int(expiration_time.timestamp())
                    }
                }
            )
            json_response = json.loads(response.text)
            ipfs_cid = json_response['IpfsHash']

            link = f"https://gateway.pinata.cloud/ipfs/{ipfs_cid}"
            return link

    except Exception as e:
        print(f"\n\nError in uploadTxt(): {e}\n\n")



def blast(fastaQuery):
    output=None
    try:

        # perform a BLAST search using NCBIWWW.qblast and generate an output file
        # containing the alignment results. The file is then uploaded to a specified location,
        # and the corresponding URL is stored in a output object.

        now = datetime.now()
        name = "BlastOutput_" + now.strftime("%m%d%y_%H_%M_%S") + "_" + str(random.randrange(0, 1000))
        outputTxtFile =  name+".txt"
        outputXmlFile = name+".xml"

        result_handle = NCBIWWW.qblast("blastp", "nr", fastaQuery.seq)

        with open("mysite/static/serve/files/"+outputXmlFile, "w") as file:
            file.write(result_handle.read())

        with open("mysite/static/serve/files/"+outputXmlFile, "r") as file:
            blast_record = NCBIXML.read(file)

        os.remove("mysite/static/serve/files/"+outputXmlFile)
        E_VALUE_THRESHOLD = 0.04

        alignmentTxt = ""
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESHOLD:
                    alignmentTxt += f"Alignment:\n"
                    alignmentTxt += f"\tSequence: {alignment.title}\n"
                    alignmentTxt += f"\tLength: {alignment.length}\n"
                    alignmentTxt += f"\tE value: {hsp.expect}\n"
                    alignmentTxt += f"\tQuery: {hsp.query[0:75]}\n"
                    alignmentTxt += f"\tMatch: {hsp.match[0:75]}\n"
                    alignmentTxt +=f"\tSbjct: {hsp.sbjct[0:75]}\n\n"

        with open("mysite/static/serve/files/"+outputTxtFile, "w") as file:
            file.write(alignmentTxt)

        return outputTxtFile

    except Exception as e:
        print("Error in blast | ",e)
        return output



class PeptideOperations():

    def __init__(self):
        self.flist={}
        self.glist={}
        self.seq=""


    ############################################
    #  Peptide Operation calculation Functions #
    ############################################

    def f1(self): #Peptide Size

        # Updates the 'Peptide Size' key in the 'flist' dictionary
        # with the length of the 'seq' string.

        self.flist.update({"Peptide Size":len(self.seq)})


    def f2(self): #Molecular Weight

        # Calculates the molecular weight of the protein sequence using the
        # 'ProteinAnalysis' class from the 'Bio.SeqUtils.ProtParam' module,
        # and updates the 'Molecular Weight' key in the 'flist' dictionary
        # with the calculated value.

        x= ProteinAnalysis(self.seq)
        ml=("%0.2f" % x.molecular_weight()+" Dalton")
        self.flist.update({"Molecular Weight":ml})


    def f3(self): #one-letter code to three-letter code

        # Converts the protein sequence from one-letter code to three-letter
        # code using the 'seq3' function from the 'Bio.SeqUtils' module,
        # and updates the 'One-to-Three Letter Code' key in the 'flist' dictionary
        # with the converted value.

        ml=seq3(self.seq, custom_map={"*": "***"},undef_code='---')
        self.flist.update({"One-to-Three Letter Code":ml})


    def f4(self): #Amino Acid Distribution

        # Calculates the frequency of each amino acid in the protein sequence,
        # and updates the 'Amino Acid Distribution' key in the 'flist' dictionary
        # with a dictionary containing the count of each amino acid in the protein
        # sequence, with keys being the full name of each amino acid.

        d={'A': 'Alanine',
           'C': 'Cysteine',
           'D': 'Aspartate',
           'E': 'Glutamate',
           'F': 'Phenylalanine',
           'G': 'Glycine',
           'H': 'Histidine',
           'I': 'Isoleucine',
           'K': 'Lysine',
           'L': 'Leucine',
           'M': 'Methionine',
           'N': 'Aspargine',
           'P': 'Proline',
           'Q': 'Glutamine',
           'R': 'Arginine',
           'S': 'Serine',
           'T': 'Threonine',
           'V': 'Valine',
           'W': 'Trypthophan',
           'Y': 'Tyrosine'}
        frequency = collections.Counter(self.seq)
        fq={}
        fq=dict(frequency)
        acid_count={}
        for i in d:
            if i in fq:
                acid_count.update({d[i]:fq[i]})
            else:
                acid_count.update({d[i]:0})
        self.flist.update({"Amino Acid Distribution":acid_count})


    def f5(self): #Amino Acid Classifier

        # Classifies the amino acids in the protein sequence into five categories
        # based on their chemical properties, and updates the 'Amino Acid Classifier'
        # key in the 'flist' dictionary with a dictionary containing the count of
        # amino acids in each category.
        #
        # The five categories are:
        #    - Aromatic (Y, W, F)
        #    - Aliphatic (A, V, I, M, L)
        #    - Basic (R, H, K)
        #    - Acidic (D, E)
        #    - Neutral (S, T, N, Q)

        frequency = collections.Counter(self.seq)
        fq={}
        fq=dict(frequency)
        Aromat=0
        Aliphat=0
        Basic=0
        Acidic=0
        Neutral=0
        cat_ac={}
        for x in fq:
            if x in ["Y","W","F"]:
                Aromat= Aromat+fq.get(x)
            elif x in ["A","V","I","M","L"]:
                Aliphat=Aliphat+fq.get(x)
            elif x in ["R","H","K"]:
                Basic=Basic+fq.get(x)
            elif x in ["D","E"]:
                Acidic=Acidic+fq.get(x)
            elif x in ["S","T","N","Q"]:
                Neutral=Neutral+fq.get(x)
        for variable in ["Aromat", "Aliphat","Basic","Acidic","Neutral"]:
            cat_ac[variable] = eval(variable)

        self.flist.update({"Amino Acid Classifier":cat_ac})


    def f6(self): #charge of peptide

        # This function calculates the charge of a peptide based on its
        # amino acid composition.

        frequency = collections.Counter(self.seq)
        fq={}
        fq=dict(frequency)
        Aromat=0
        Aliphat=0
        Basic=0
        Acidic=0
        Neutral=0
        cat_ac={}
        for x in fq:
            if x in ["Y","W","F"]:
                Aromat= Aromat+fq.get(x)
            elif x in ["A","V","I","M","L"]:
                Aliphat=Aliphat+fq.get(x)
            elif x in ["R","H","K"]:
                Basic=Basic+fq.get(x)
            elif x in ["D","E"]:
                Acidic=Acidic+fq.get(x)
            elif x in ["S","T","N","Q"]:
                Neutral=Neutral+fq.get(x)
        charge= (Basic-Acidic)
        if charge>1:
            ml=f"Charge of the peptide is {charge}."
        elif charge ==0:
            ml="Peptide is neutral"
        else:
            ml=f"Charge of the peptide is {charge}."

        self.flist.update({"Peptide Charge":ml})


    def f7(self): #aromaticity

        # Calculates the percentage of aromatic amino acids in a protein
        # sequence, and updates the 'flist' dictionary with a key-value
        # pair containing the result.

        x=ProteinAnalysis(self.seq)
        ml=("%0.2f" % x.aromaticity())
        self.flist.update({"The percentage Aromaticity of the peptide is":ml})


    def f8(self): #half life and gravy index

        # Predicts the half-life of a protein sequence by calculating
        # its instability index, and updates the 'flist' dictionary
        # with a key-value pair containing the prediction.

        x = ProteinAnalysis(self.seq)
        ml = float("%0.2f" % x.instability_index())
        gravy = float("%0.2f" % x.gravy())

        if ml > 40:
            self.flist.update({"Half Life and Gravy Index": f"Half life of peptide is greater than 40 therefore peptide is unstable. Gravy index value is {gravy}"})
        else:
            self.flist.update({"Half Life and Gravy Index": f"Half life of peptide is smaller than 40 therefore peptide is stable. Gravy index value is {gravy}"})


    def f9(self): #secondary structure

        # Predicts the secondary structure of a protein sequence by
        # classifying each amino acid in the sequence as either alpha-helix,
        # beta-sheet, turn, or irregular structure, and updates the 'flist'
        # dictionary with a key-value pair containing the prediction.

        ml=[]
        code=[]
        for x in self.seq:
            code.append(x)
        for i in code:
            if i in ['A','C','L','M','Q','E','H','K']:
                ml.append("ALPHA")
            elif i in ['V','I','F','Y','W','T','R']:
                ml.append("BETA")
            elif i in ['G','S','P','D','N']:
                ml.append("TURN")
            else:
                ml.append("IS")
        self.flist.update({"Predicted Secondary Structure":ml})


    def f10(self): #Molar Extinction Coefficient

        # Calculates the molar extinction coefficient of the protein sequence,
        # specifically for the reduced form of the protein
        # (i.e., with no disulfide bonds between cysteine residues).
        # The result is added to the 'Molar Extinction Coefficient' key in the 'flist' dictionary.

        x=ProteinAnalysis(self.seq)
        epsilon_prot=x.molar_extinction_coefficient()  # [reduced, oxidized]
        ml= f"{epsilon_prot[0]} with reduced number of cysteines"
        self.flist.update({"Molar Extinction Coefficient":ml})


    def f11(self): #Disulphide Bridges Calculation

        # Calculates the number of disulfide bonds present in the protein sequence,
        # specifically for the oxidized form of the protein
        # (i.e., with disulfide bonds between cysteine residues).
        # The result is added to the 'Disulphide Bridges Calculation' key in the 'flist' dictionary.

        x=ProteinAnalysis(self.seq)
        epsilon_prot=x.molar_extinction_coefficient()  # [reduced, oxidized]
        ml=f"The number of disulphide bridges present is {epsilon_prot[1]}"
        self.flist.update({"Disulphide Bridges":ml})


    def f12(self): #Binding Potential Calculator

        # Calculates the Boman index value of the given protein sequence using the
        # Global Descriptor tool and classifies the peptide based on its binding potential.
        # Updates the 'Binding Potential Calculator' key in the 'flist' dictionary with a string
        # containing the Boman index value and the binding potential classification.

        desc=GlobalDescriptor(self.seq)
        bomanvalue=desc.boman_index()
        val=desc.descriptor
        value = ''
        for i in val:
            for j in i:
                value = value + str(j) + ', '
        if val > 2.48:
            ml=f"Boman index value of peptide is {value}It has high binding potential"
        else:
            ml=f"Boman index value of peptide is {value}Low binding potential"
        self.flist.update({"Binding Potential":ml})


    def f17(self): #isoelectricpoint

        pi_7 = float("%0.2f" % IsoelectricPoint(self.seq).pi())
        charge_at_pH_7 = float("%0.2f" % IsoelectricPoint(self.seq).charge_at_pH(7.0))
        charge_at_pH_4_53 = float("%0.2f" % ProteinAnalysis(self.seq).charge_at_pH(4.53))
        charge_at_pH_10_53 = float("%0.2f" % ProteinAnalysis(self.seq).charge_at_pH(10.53))

        self.flist.update({"Isoelectric Point and Charge":
                        f"Isoelectric point of peptide is {pi_7}. "
                        f"Charge at pH 7 is {charge_at_pH_7}. "
                        f"Charge at pH 4.53 is {charge_at_pH_4_53}. "
                        f"Charge at pH 10.53 is {charge_at_pH_10_53}."})




    ###################################################
    #  Graphical Representations of Peptide Functions #
    ###################################################

    def f13(self): #hydropathy plot

        # Generates a hydropathy plot for the given peptide sequence using the
        # Kyte-Doolittle hydrophobicity scale. Returns the URL of the uploaded image.

        try:
            kd = {  'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
                    'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
                    'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
                    'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2
                 }
            values = []
            num_residues=len(self.seq)
            for residue in self.seq:
                values.append(kd[residue])
            x_data = range(1, num_residues+1)
            plt.plot(x_data, values, linewidth=1.0)
            plt.axis(xmin = 1, xmax = num_residues)
            plt.xlabel("Residue Number")
            plt.ylabel("Hydrophobicity")

            # Generate a unique file name for the plot
            now=datetime.now()
            file_name = "hydropathy" + now.strftime("%m%d%y_%H_%M_%S") + "_" + str(random.randrange(0, 1000)) + ".png"
            #plt.gcf().set_size_inches(10, 5)
            plt.savefig("mysite/static/serve/images/"+file_name, bbox_inches='tight')
            plt.clf()


            self.glist.update({"Hydropathy Plot": file_name})
            #os.remove(file_name)

        except Exception as e:
            print(f"\n\nError in f13(): {e}\n\n")


    def f14(self): #Helical Wheel Projection

        # Generates a helical wheel projection for the given peptide sequence using the 'helical_wheel'
        # function from the 'modlamp' library, and uploads the generated image to imgbb.com.
        # The function updates the 'Helical Wheel Projection for Peptide' key in the 'glist' dictionary
        # with the URL of the uploaded image, and deletes the local image file after the upload is complete.

        now=datetime.now()
        file_name = "helical_Wheel_" + now.strftime("%m%d%y_%H_%M_%S") + "_" + str(random.randrange(0, 1000)) + ".png"

        try:
            # Generate the helical wheel projection and save it as a PNG file
            helical_wheel(self.seq, colorcoding='amphipathic', lineweights=True, filename="mysite/static/serve/images/"+file_name, seq=False, moment=True)


            self.glist.update({"Helical Wheel Projection": file_name})
            #os.remove(file_name)

        except Exception as e:
            print(f"\n\nError in f14(): {e}\n\n")


    def f15(self): # Amino acid distribution plot
        # Amino acid distribution plot

        try:
            plt.clf()

            d = {'A': 'Alanine',
                'C': 'Cysteine',
                'D': 'Aspartate',
                'E': 'Glutamate',
                'F': 'Phenylalanine',
                'G': 'Glycine',
                'H': 'Histidine',
                'I': 'Isoleucine',
                'K': 'Lysine',
                'L': 'Leucine',
                'M': 'Methionine',
                'N': 'Aspargine',
                'P': 'Proline',
                'Q': 'Glutamine',
                'R': 'Arginine',
                'S': 'Serine',
                'T': 'Threonine',
                'V': 'Valine',
                'W': 'Trypthophan',
                'Y': 'Tyrosine'}

            frequency = collections.Counter(self.seq)
            fq={}
            fq=dict(frequency)
            acid_count={}
            for i in d:
                if i in fq:
                    acid_count.update({d[i]:fq[i]})

            print(acid_count)

            plt.bar(acid_count.keys(), acid_count.values())
            plt.xlabel('Amino acids')
            plt.ylabel('Amino acid frequency')
            plt.gca().yaxis.set_major_locator(ticker.MaxNLocator(integer=True)) # only show integer y-axis ticks
            if len(acid_count.values())>10:
                plt.xticks(rotation=90) # rotate x-axis tick labels by 90 degrees
            #plt.show()

            # Generate a unique file name for the plot
            now=datetime.now()
            file_name = "Distribution_" + now.strftime("%m%d%y_%H_%M_%S") + "_" + str(random.randrange(0, 100000)) + ".png"
            plt.savefig("mysite/static/serve/images/"+file_name, bbox_inches='tight')
            plt.clf()

            self.glist.update({"Distribution of Amino Acid": file_name})

        except Exception as e:
            print(f"\n\nError in f15(): {e}\n\n")


    def f16(self): #Amino Acid classification

        # Categorise amino acids and plotting pie chart and predicting overall charge on peptide
        try:
            Aromat = self.seq.count("Y" or "y")+ self.seq.count("W" or "w") + self.seq.count("F" or "f")
            Aliphat = self.seq.count("A" or "a") + self.seq.count("V" or "v")+ self.seq.count("I" or "i") + self.seq.count("M" or "m")+ self.seq.count("L" or "l")
            Basic = self.seq.count("R" or "r")+ self.seq.count("H" or "h") + self.seq.count("K" or "k")
            Acidic = self.seq.count("D" or "d") + self.seq.count("E" or "e")
            Neutral = self.seq.count("S" or "s") +  self.seq.count("T" or "t") + self.seq.count("N" or "n") + self.seq.count("Q" or "q")

            y = np.array([Aliphat, Basic, Acidic, Neutral, Aromat])
            mylabels = ["Aliphat", "Basic", "Acidic", "Neutral", "Aromat"]

            plt.pie(y, labels = mylabels,autopct='%1.1f%%')
            plt.legend(title = "Amino acid Classification:")

            # Generate a unique file name for the plot
            now = datetime.now()
            file_name = "classification_" + now.strftime("%m%d%y_%H_%M_%S") + "_" + str(random.randrange(0, 1000)) + ".png"

            plt.savefig("mysite/static/serve/images/"+file_name, bbox_inches='tight')
            plt.clf()

            self.glist.update({"Amino Acid classification": file_name})


        except Exception as e:
            print(f"\n\nError in f16(): {e}\n\n")


    def check(self, seq:str, check_list:list) -> dict:

        self.seq = seq.upper()
        result = []

        for i in check_list:
            if i== 'f1':
                self.f1()
            elif i== 'f2':
                self.f2()
            elif i== 'f3':
                self.f3()
            elif i== 'f4':
                self.f4()
            elif i== 'f5':
                self.f5()
            elif i== 'f6':
                self.f6()
            elif i== 'f7':
                self.f7()
            elif i== 'f8':
                self.f8()
            elif i== 'f9':
                self.f9()
            elif i== 'f10':
                self.f10()
            elif i== 'f11':
                self.f11()
            elif i== 'f12':
                self.f12()
            elif i== 'f13':
                self.f13()
            elif i== 'f14':
                self.f14()
            elif i== 'f15':
                self.f15()
            elif i== 'f16':
                self.f16()
            elif i== 'f17':
                self.f17()

        result.append(self.flist)
        result.append(self.glist)
        self.flist={}
        self.glist={}
        self.seq=""
        return(result)
