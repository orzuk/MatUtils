function indep = get_dag_indep( G)
% Get ALL the independecies from a graph G 
% The number we check is very high : 4^n/2, so use this
% only for small values of n. The coding is as follows : 
% X is indep of Y given Z wrt DAG G?
%
% X : 0,  Y : 1,  Z : 2,  V \ {X U Y U Z} : 3
%
% Note: This function is not needed so much, since it is enough to check
% the Local Markov properties!!! 
%

n = size(G,1);

indep = zeros(1,4^n); % Here we duplicate everything twice - Lo Nora.

for i=0:4^n-1

    digs = mod(floor(i ./ 4.^[0:n-1]), 4);
    X = find(digs==0); Y = find(digs==1); Z = find(digs==2);
    
    if( (~isempty(X)) && (~isempty(Y)) )
        if(X(1) < Y(1))
            indep(i+1) = dsep(X, Y, Z, G);
        end
    end

    % % % % % %     % Get X,Y and Z
    % % % % % %    XYZ={}; XYZ{1}=[]; XYZ{2}=[]; XYZ{3}=[]; XYZ{4}=[];
    % % % % % %
    % % % % % %    for j=0:n-1 % Get the digits
    % % % % % %        dig = mod(floor(i / 4^j), 4);
    % % % % % %
    % % % % % %        XYZ{dig+1} = [XYZ{dig+1} j+1];
    % % % % % %    end
    % % % % % %
    % % % % % % %    if(i < 11)
    % % % % % % %        i
    % % % % % % %        X= XYZ{1}
    % % % % % % %        Y= XYZ{2}
    % % % % % % %        Z= XYZ{3}
    % % % % % % %        ind_i  = indep(i+1)
    % % % % % % %        G
    % % % % % % %    end
    % % % % % %
    % % % % % %    indep(i+1) = dsep(XYZ{1}, XYZ{2}, XYZ{3}, G);

end
