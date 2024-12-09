import plotly.express as px
import pandas as pd

if __name__ == "__main__":
    #Hyper Fold Enrichment, Hyper Total Genes/Regions
    data = {"GO_CC_pre-autophagosomal_structure_membrane": (15.3910, 16),
            "GO_BP_myo-inositol_transport": (33.3333, 6),
            "GO_CC_tRNA_(m1A)_methyltransferase_complex": (72.0370, 2),
            "GO_CC_A_band": (9.4595, 37), 
            "GO_CC_M_band": (12, 25)}
    

    pathways = [key for key in data]
    Hyper_FE = [data[key][0] for key in data]
    Hyper_regions = [data[key][1] for key in data]

    df = pd.DataFrame(dict(Pathway=pathways, Fold_Enrichment=Hyper_FE, Genes=Hyper_regions))
    
    fig = px.scatter(df, x="Fold_Enrichment", y="Pathway", size="Genes",
                 title="Pathway Fold Enrichment",
                 labels={"salary":"Annual Salary (in thousands)"} # customize axis label
                )
    fig.show()