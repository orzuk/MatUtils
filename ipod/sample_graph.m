% Sample a graph, which means taking an existing
% graph and add some noise on it, for example miss some
% of the edges/node or explore the graph partially etc.
function [G_samp]=sample_graph(G, sampling_type, param1, param2)
SAMP_EDGES = 0; SAMP_NODES = 1;  % sampling flags 

N = length(G);

% Simple edge sampling: each edge which is ON is ruined with prob. p_on,
% and each edge which is OFF, is turned on with prob. p_off
if(sampling_type == SAMP_EDGES)
    p_off = param1; p_on = param2;
%     
% % % %     q = sqrt(p_off); % Compensate for randomizing twice
% % % %    
% % % %     E_s=rand(N) < q;
% % % %     E_s=E_s.*E_s';
% % % %     for i=1:N
% % % %         E_s(i,i)=0;
% % % %     end
% % % %     E_s = sparse(E_s); % turn into sparse

   
    E_off = generate_random_graph(0, N, p_off); % Generate the edges which should be off    
    G_samp = max(G-E_off,0);
    
    E_on = generate_random_graph(0, N, p_on); % Generate the edges which should be off
    G_samp = or(G_samp, E_on); % Slight inaccuracy here. If p_off turned it off then p_on can still turn it on
    
    
    
 %   E_s = sprandsym(N,p); 
    
% % % % %     if(remove_only)
% % % % %         G_samp = and(xor(E_s, G),G); 
% % % % %     else
% % % % %         G_samp = xor(E_s, G); 
% % % % %     end
    
    % Make diagonal zero
    for i=1:N 
        G_samp(i,i)=0;
    end
    
    return;
    
end

if(sampling_type == SAMP_NODES)
   p = param; 
   
   ns = find(rand(1,N) < 1-p); % We throw a fraction p of the nodes   
   G_samp = G(ns,ns);
    
   return;
end
    
    

