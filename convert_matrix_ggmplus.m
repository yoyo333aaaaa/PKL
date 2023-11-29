N=length(X); 
X_New=[]; 
for i=1:N; 
    X_New=[X_New;X(:,i)]; 
end 
X_New 

M=length(Y); 
Y_New=[]; 
for i=1:N; 
    Y_New=[Y_New;Y(:,i)]; 
end 
Y_New 

O=length(Z); 
Z_New=[]; 
for i=1:N; 
    Z_New=[Z_New;Z(:,i)]; 
end 
Z_New