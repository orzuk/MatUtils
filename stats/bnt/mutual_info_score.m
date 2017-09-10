function score = mutual_info_score(i,si,j,sj,data)
% G = mutual_info_score(i,si,j,sj,data)
% Only for tabular node which values are 1,2,...,size .
% si is size of node i, sj is size of node j.
% data(i,m) is the node i in the case m.
% 
% cf C.CHOW & C.LIU : Approximating discrette probability distributions with dependences trees.
%
% olivier.francois@insa-rouen.fr, philippe.leray@insa-rouen.fr

[n N]=size(data);
Nj=hist(data(j,:),1:sj);
Ni=hist(data(i,:),1:si);
NiNj=Ni'*Nj;

for k=1:si
 ind=find(data(i,:)==k) ;
 Nij(k,:) = hist(data(j,ind),1:sj);
end

% sommons les valeurs non-infinies:
ind=find(NiNj~=0 & Nij~=0);
score=sum(sum(Nij(ind).*log(Nij(ind)./NiNj(ind))/N));

% ajoutons un à priori de Dirichlet d'exposant (1...1) pour eviter
% d'avoir des log(0)=-inf, le score devient alors :
% score=sum(sum(((Nij+1)./(N+si*sj)).*log(((Nij+1)*(N+si)*(N+sj))./(si*sj.*NiNj))));

