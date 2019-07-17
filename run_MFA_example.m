File1 = {'sample input file 1.xlsx'};
File2 = {'sample input file 2.xlsx'};
File3 = {'sample input file 3.xlsx'};

[emumodel]=createmumodel(File1, File2, File3);

[res]=flxestimate(emumodel);

[res,impres]=confintestimate(res,emumodel);