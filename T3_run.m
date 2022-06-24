%% Preparación, Preguntas 1 2 3
obj = T3_obj('docvis.xlsx');

%% Pregunta 4
beta_est_NewtonR = obj.MetNum("NR",100,[0 0 0 0 0],obj.regresores,obj.Y);
%% Pregunta 5
beta_est = obj.MetNum("BHHH",1000,[0 0 0 0 0],obj.regresores,obj.Y);

%% Pregunta 6
beta_est_rest = obj.NR_R(100,obj.Y,0);

%% Pregunta 7 Test Wald
Wald_st = obj.TestWald(obj.regresores,obj.Y,beta_est_NewtonR);
%% Pregunta 7 Test Multiplicadores de Lagrange
lm_st =obj.MultLag(obj.regresores,obj.Y,beta_est_rest);
%% Pregunta 7 Test de Likelihood Ratio
lkrat_st = obj.Lratio(obj.regresores,obj.Y,beta_est_NewtonR,beta_est_rest);
