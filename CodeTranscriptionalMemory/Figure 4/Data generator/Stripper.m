ClonekeyN = {'Rasa12_S1_STRSa1','Pgk_S1_STRSa1','Rbpj_S1_STRSa1','Jarid2_S1_STRSa1','Jam2_S1_STRSa1','Nono_S1_STRSa1','Tfpd2_S1_STRSa1','Dstn_S1_STRSa1','Sprouty_S1_STRSa1'};

for z3 = 1:length(ClonekeyN)
    load([char(ClonekeyN(z3))])
    ClonekeyN = {'Rasa12_S1_STRSa1','Pgk_S1_STRSa1','Rbpj_S1_STRSa1','Jarid2_S1_STRSa1','Jam2_S1_STRSa1','Nono_S1_STRSa1','Tfpd2_S1_STRSa1','Dstn_S1_STRSa1','Sprouty_S1_STRSa1'};
    save([char(ClonekeyN(z3))],'X','x','yTOT','gamma','R0R0coeff')
  
end