"""Publication-facing identifiers for the 15 biological RNA-seq samples."""

PUBLIC_SAMPLES = {
    "TGN1003AW1": "Cheng_TGN_L1_W1",
    "TGN1003AW2": "Cheng_TGN_L1_W2",
    "TGN1003AW3": "Cheng_TGN_L1_W3",
    "TGN1003AW4": "Cheng_TGN_L1_W4",
    "TGN251028BW1": "Cheng_TGN_Pilot_W1",
    "TGN1028BW1": "Cheng_TGN_L2_W1",
    "TGN1028BW2": "Cheng_TGN_L2_W2",
    "TGN1028BW3": "Cheng_TGN_L2_W3",
    "TGN1028BW4": "Cheng_TGN_L2_W4",
    "TGN1215AvXW1": "Cheng_TGN_L3_W1",
    "TGN1215AvXW2": "Cheng_TGN_L3_W2",
    "TGN1215AvXW3": "Cheng_TGN_L3_W3",
    "TGN1215AvXW4": "Cheng_TGN_L3_W4",
    "hDRG": "Cheng_hDRG1",
    "TAK2204464A": "Cheng_hDRG2",
}

def public_name(source_sample):
    source_sample = source_sample.removeprefix("C1_").removeprefix("C2_")
    return PUBLIC_SAMPLES[source_sample]
