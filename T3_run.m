%% Preparación, Preguntas 1 2 3
obj = T3_obj('docvis.xlsx');

%% Pregunta 4
beta_est = obj.MetNum("NR",100,[0 0 0 0 0],obj.regresores,obj.Y);
%% Pregunta 5
beta_est = obj.MetNum("BHHH",1000,[0 0 0 0 0],obj.regresores,obj.Y);

%% Pregunta 6
beta_est = obj.NR_R(100,obj.Y,0);

