
%%author: Marius Leordeanu
% modified: May 16, 2011

%Utility: example script for running IPFP stand-alone and IPFP
%         combined with spectral matching on the well-known "tho30" QAP problem 


%% A and B obtained from tho30.dat


A = [
  0   1   2   3   4   5   6   7   8   9   1   2   3   4   5   6   7   8   9  10   2   3   4   5   6   7   8   9  10  11; 
  1   0   1   2   3   4   5   6   7   8   2   1   2   3   4   5   6   7   8   9   3   2   3   4   5   6   7   8   9  10;
  2   1   0   1   2   3   4   5   6   7   3   2   1   2   3   4   5   6   7   8   4   3   2   3   4   5   6   7   8   9; 
  3   2   1   0   1   2   3   4   5   6   4   3   2   1   2   3   4   5   6   7   5   4   3   2   3   4   5   6   7   8 ;
  4   3   2   1   0   1   2   3   4   5   5   4   3   2   1   2   3   4   5   6   6   5   4   3   2   3   4   5   6   7 ;
  5   4   3   2   1   0   1   2   3   4   6   5   4   3   2   1   2   3   4   5   7   6   5   4   3   2   3   4   5   6 ;
  6   5   4   3   2   1   0   1   2   3   7   6   5   4   3   2   1   2   3   4   8   7   6   5   4   3   2   3   4   5 ;
  7   6   5   4   3   2   1   0   1   2   8   7   6   5   4   3   2   1   2   3   9   8   7   6   5   4   3   2   3   4 ;
  8   7   6   5   4   3   2   1   0   1   9   8   7   6   5   4   3   2   1   2  10   9   8   7   6   5   4   3   2   3 ;
  9   8   7   6   5   4   3   2   1   0  10   9   8   7   6   5   4   3   2   1  11  10   9   8   7   6   5   4   3   2 ;
  1   2   3   4   5   6   7   8   9  10   0   1   2   3   4   5   6   7   8   9   1   2   3   4   5   6   7   8   9  10 ;
  2   1   2   3   4   5   6   7   8   9   1   0   1   2   3   4   5   6   7   8   2   1   2   3   4   5   6   7   8   9 ;
  3   2   1   2   3   4   5   6   7   8   2   1   0   1   2   3   4   5   6   7   3   2   1   2   3   4   5   6   7   8 ;
  4   3   2   1   2   3   4   5   6   7   3   2   1   0   1   2   3   4   5   6   4   3   2   1   2   3   4   5   6   7 ;
  5   4   3   2   1   2   3   4   5   6   4   3   2   1   0   1   2   3   4   5   5   4   3   2   1   2   3   4   5   6 ;
  6   5   4   3   2   1   2   3   4   5   5   4   3   2   1   0   1   2   3   4   6   5   4   3   2   1   2   3   4   5 ;
  7   6   5   4   3   2   1   2   3   4   6   5   4   3   2   1   0   1   2   3   7   6   5   4   3   2   1   2   3   4 ;
  8   7   6   5   4   3   2   1   2   3   7   6   5   4   3   2   1   0   1   2   8   7   6   5   4   3   2   1   2   3 ;
  9   8   7   6   5   4   3   2   1   2   8   7   6   5   4   3   2   1   0   1   9   8   7   6   5   4   3   2   1   2 ;
 10   9   8   7   6   5   4   3   2   1   9   8   7   6   5   4   3   2   1   0  10   9   8   7   6   5   4   3   2   1 ;
  2   3   4   5   6   7   8   9  10  11   1   2   3   4   5   6   7   8   9  10   0   1   2   3   4   5   6   7   8   9 ;
  3   2   3   4   5   6   7   8   9  10   2   1   2   3   4   5   6   7   8   9   1   0   1   2   3   4   5   6   7   8 ;
  4   3   2   3   4   5   6   7   8   9   3   2   1   2   3   4   5   6   7   8   2   1   0   1   2   3   4   5   6   7 ;
  5   4   3   2   3   4   5   6   7   8   4   3   2   1   2   3   4   5   6   7   3   2   1   0   1   2   3   4   5   6 ;
  6   5   4   3   2   3   4   5   6   7   5   4   3   2   1   2   3   4   5   6   4   3   2   1   0   1   2   3   4   5 ;
  7   6   5   4   3   2   3   4   5   6   6   5   4   3   2   1   2   3   4   5   5   4   3   2   1   0   1   2   3   4 ;
  8   7   6   5   4   3   2   3   4   5   7   6   5   4   3   2   1   2   3   4   6   5   4   3   2   1   0   1   2   3 ;
  9   8   7   6   5   4   3   2   3   4   8   7   6   5   4   3   2   1   2   3   7   6   5   4   3   2   1   0   1   2 ;
 10   9   8   7   6   5   4   3   2   3   9   8   7   6   5   4   3   2   1   2   8   7   6   5   4   3   2   1   0   1 ;
 11  10   9   8   7   6   5   4   3   2  10   9   8   7   6   5   4   3   2   1   9   8   7   6   5   4   3   2   1   0 ;
];
 
 
B = [ 
  0 126   0   3 224   0   0 241  21 175 133   0   0 162   0   0   0   0 155   0   0   0 114 243   0  44  48  94  43 225; 
126   0   0   0   0 137  14   0   0  13   0 245   0 128  16   0 239   0 108   0  28   0   0   0 111   0 232 144  35  44 ;
  0   0   0   0 242 180 151   0   0   0 183   0 127  26  10   0  99 131 155   6   0   0   0  29   0   0   0   0   0  45 ;
  3   0   0   0 231  15  17   0 159   0   0 125 138  44   0 210 133   0   0   0   0   0   0 150   0 207 232   0   0   0 ;
224   0 242 231   0   0 196  71   0   0  64   0  26   0   0 198  94   0   0   0   0   0   0 176   0   0  31   0 105 114 ;
  0 137 180  15   0   0   0   0 204 169   0 247 195  96 121   0   0 203   0  68   0   0   0   0  51 214   0  23  24   0 ;
  0  14 151  17 196   0   0   0   0   0 188   0   0 210   0 132  11  59   0   0  37 238   0 150   0 136 108   0   0   0 ;
241   0   0   0  71   0   0   0   0   0   0 156   0 178   0 143 208   0   0   0 115  73   0   0 167  91   0 209   0 111 ;
 21   0   0 159   0 204   0   0   0 206   0 149   0   0   0   0   0   0   0   0   0   0   0   0 201   0 210  36   4   0 ;
175  13   0   0   0 169   0   0 206   0   5   0 127   0   0   0   0  40   0  20 218   0 112   0   0 164 146  50   0 236 ;
133   0 183   0  64   0 188   0   0   5   0   0  12 120  74   0   0  74  25  58  86   0   0 190   0  81 162   3   0   0 ;
  0 245   0 125   0 247   0 156 149   0   0   0   0   1   0   0   0   0   0  39   0 110 151   0  68 197   0   0   0  89 ;
  0   0 127 138  26 195   0   0   0 127  12   0   0   0  96  52   0 182   1   0 104  82 146  64 189  17 231   0   0   0 ;
162 128  26  44   0  96 210 178   0   0 120   1   0   0   0   0   0   0 123  84 127 198 159   8   0  61   0  61   0   0 ;
  0  16  10   0   0 121   0   0   0   0  74   0  96   0   0   0  93 216  44  12   0 229   0   0   0 141  21 114   0 157 ;
  0   0   0 210 198   0 132 143   0   0   0   0  52   0   0   0  98  13   0  69 242  22 193   0   0  36  16  80  47   0 ;
  0 239  99 133  94   0  11 208   0   0   0   0   0   0  93  98   0   0   0   0  93  95  81   0   0 234 126 170  23  40 ;
  0   0 131   0   0 203  59   0   0  40  74   0 182   0 216  13   0   0   0   0 119   0   0  18  79   0  17  49   0   0 ;
155 108 155   0   0   0   0   0   0   0  25   0   1 123  44   0   0   0   0   0   0   0   0 139   0 147   0  28 133  82 ;
  0   0   6   0   0  68   0   0   0  20  58  39   0  84  12  69   0   0   0   0   0 236  86   0   0   0   0 172   0   7 ;
  0  28   0   0   0   0  37 115   0 218  86   0 104 127   0 242  93 119   0   0   0 183 214   0   0   0 100   0   0  60 ;
  0   0   0   0   0   0 238  73   0   0   0 110  82 198 229  22  95   0   0 236 183   0  75 113 209 211   0  87   0  61 ;
114   0   0   0   0   0   0   0   0 112   0 151 146 159   0 193  81   0   0  86 214  75   0 140   0  49   0  44   7   0 ;
243   0  29 150 176   0 150   0   0   0 190   0  64   8   0   0   0  18 139   0   0 113 140   0 203 232 214 121   0   0 ;
  0 111   0   0   0  51   0 167 201   0   0  68 189   0   0   0   0  79   0   0   0 209   0 203   0   0 153 200   0   0 ;
 44   0   0 207   0 214 136  91   0 164  81 197  17  61 141  36 234   0 147   0   0 211  49 232   0   0   0 115   0 103 ;
 48 232   0 232  31   0 108   0 210 146 162   0 231   0  21  16 126  17   0   0 100   0   0 214 153   0   0  62 159   0 ;
 94 144   0   0   0  23   0 209  36  50   3   0   0  61 114  80 170  49  28 172   0  87  44 121 200 115  62   0 229  90 ;
 43  35   0   0 105  24   0   0   4   0   0   0   0   0   0  47  23   0 133   0   0   0   7   0   0   0 159 229   0   0 ;
225  44  45   0 114   0   0 111   0 236   0  89   0   0 157   0  40   0  82   7  60  61   0   0   0 103   0  90   0   0 ;  


];


nNodes = 30;
nLabels = 30;

M = zeros(nNodes*nLabels);

labels = zeros(1, nNodes*nLabels);
nodes = zeros(1, nNodes*nLabels);


el = 0;

for node =  1:nNodes
    for label = 1:nLabels
   
       el = el + 1;
       
       nodes(el)  = node;
       labels(el) = label;
        
    end
end


for i = 1:length(nodes)
    for j = 1:length(labels)
   
        M(i,j) =  A(nodes(i), nodes(j))*B(labels(i), labels(j));
        
    end
end


sol0 = ones(length(nodes),1);

sol0 = sol0/norm(sol0);


%% test ipfp stand-alone starting from a flat/uniform solution

[sol_ipfp, stats_ipfp]  = ipfp_gm(M, sol0, labels, nodes);

score_ipfp = sol_ipfp'*M*sol_ipfp/2;


%% test ipfp starting from the solution returned by spectral matching

[sol_ipfp_sm, stats_ipfp_sm]  = spectral_matching_ipfp(M, labels, nodes);

score_ipfp_sm = sol_ipfp_sm'*M*sol_ipfp_sm/2;


%% test the second ipfp matlab function stand alone

D = zeros(length(sol0), 1);

[sol_ipfp2, x_opt, score, score_sol]  = ipfp(M, D, sol0, labels, nodes, 50);

score_ipfp2 = sol_ipfp2'*M*sol_ipfp2/2;





return;









