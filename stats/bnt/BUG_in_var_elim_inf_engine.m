% Small test for detecting a bug in BNT var_elim_inf_engine !!!! 
b = [];
dag = [0 0; 1 0];
b = mk_bnet(dag, [2 2]); % make a two-by-two net

b.CPD{1} = tabular_CPD(b, 1, 'CPT', [0.3 0.7; 0.8 0.2]);
b.CPD{2} = tabular_CPD(b, 2, 'CPT', [0.1 0.9]);

P = zuk_bnet_to_probs(b)  %Get joint probability distribution
evidence = cell(1, 2);
engine = var_elim_inf_engine(b);
engine = enter_evidence(engine, evidence);
m1 = marginal_nodes(engine,1); 
m2 = marginal_nodes(engine,2);
m1.T
m2.T


