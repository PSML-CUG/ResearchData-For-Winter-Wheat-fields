%   ResearchData_Package Name：Temporal variation of scattering powers from Freeman-Durden Decomposition  over Winter Wheat fields
%  Author: Wenxin Xue (wenxin@cug.edu.cn), China University of Geosciences (Wuhan)
%  Tutor: Qinghua Xie (xieqh@cug.edu.cn), China University of Geosciences (Wuhan)
%  Date: March 2025
% Copyright (c) 2025 by Wenxin Xue and Qinghua Xie, China University of Geosciences (Wuhan)

%In this study, C-band RADARSAT-2 imagery was utilized. (Please note that this data is not subject to copyright restrictions; therefore, I will provide the power values of its components following Freeman-Durden decomposition.) 

%A T3 matrix generated from the May 9, 2019 image ('T3_0509.mat') is provided for code testing.

%Step1: Involved preprocessing the data in SNAP 9.0 and generating the coherence matrix (T3). 

%Step 2: The T3 matrix was imported into MATLAB, where the double-bounce, volume scattering power, and surface scattering power for each image were computed using a combination of Freeman-Durden three-component decomposition and non-negative eigenvalue decomposition (refer to Section 4.1 of the paper).

%Step 3: These results were integrated into a three-dimensional matrix (line*column*3), yielding the file 'freeman_pspdpv_date.mat'.

%%%   The programs contained in this package are granted free of charge for
   research and education purposes only. Scientific results produced using
   the software provided shall acknowledge the use of this implementation
   provided by us. If you plan to use it for non-scientific purposes,
   don't hesitate to contact us. Because the programs are licensed free of
   charge, there is no warranty for the program, to the extent permitted
   by applicable law. except when otherwise stated in writing the
   copyright holders and/or other parties provide the program "as is"
   without warranty of any kind, either expressed or implied, including,
   but not limited to, the implied warranties of merchantability and
   fitness for a particular purpose. the entire risk as to the quality and
   performance of the program is with you. should the program prove
   defective, you assume the cost of all necessary servicing, repair or
   correction. In no event unless required by applicable law or agreed to
   in writing will any copyright holder, or any other party who may modify
   and/or redistribute the program, be liable to you for damages,
   including any general, special, incidental or consequential damages
   arising out of the use or inability to use the program (including but
   not limited to loss of data or data being rendered inaccurate or losses
   sustained by you or third parties or a failure of the program to
   operate with any other programs), even if such holder or other party
   has been advised of the possibility of such damages.
