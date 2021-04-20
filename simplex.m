format short
clear all 



RUN=true;
choice=3;
while (choice ~= 1 && choice ~=2) 
choice=input("choisir 1 pour maximisation 2 pour minimisation");
end
if (choice ==1)
NbreVariables=input('Donner le nombre des variables');
disp('donner les coefficients du Max dans l ordre');
for i=1:NbreVariables
    C(i)=input(strcat('C [', num2str(i),']=   '));
end
Ineq=input('Donner le nombre des inequations');
disp ('Donner les coefficients des varialbles dans les inequations dans l ordre');
for i=1:Ineq
    for j=1:NbreVariables
            Info(i,j)=input(strcat('Info [', num2str(i),',',num2str(j),']=   '));
    end
end
disp('donner les coefficients à droite de signe des inequations dans l ordre');
for i=1:Ineq
    b(i,1)=input(strcat('Coeff [', num2str(i),']=   '));
end
    
   s=eye(size(Info,1));

A=[Info s b];
K=["Val";"Zj-Cj"];
for i=1:Ineq
    K(i+2,1)=strcat("e_",num2str(i));
end
for i=1:NbreVariables
    J(i)=strcat("x_",num2str(i));
end
for i=NbreVariables+1:Ineq+NbreVariables
    J(i)=strcat("e_",num2str(i-NbreVariables));
end
J(NbreVariables+Ineq+1)="bj";

K1=array2table(K);

Cost=zeros(1,size(A,2));
Cost(1:NbreVariables)=C;

BV= NbreVariables+1:1:size(A,2)-1;

ZjCj=Cost-Cost(BV)*A;
ZCj=[ZjCj ; A];
B=[J ; ZCj];
MySimplexTable=array2table(B);
MySimplexTable=[K1 MySimplexTable];
disp(MySimplexTable);
    
while RUN
if any(ZjCj>0);
     ZC=ZjCj(1:end-1);
     [EnterCol,pvt_col]=max(ZC);
     fprintf ('Ke min de Zj-Cj est %d de la colonne %d \n',EnterCol,pvt_col);
     fprintf('Variable entrante est '); disp(MySimplexTable(1,pvt_col+1));
    sol= A(:,end);
    Column=A(:,pvt_col);
    if all(A(:,pvt_col)<=0)
        error('Impossible de resoudre, toutes les valeurs entrees <=0 dans la colonne %d',pvt_col);
    else
    for i=1:size(Column,1)
        if Column(i)>0
            ratio(i)=sol(i)./Column(i);
        else
            ratio(i)=inf;
        end
    end
    [MinRatio,pvt_row]=min(ratio);
    disp('Variable sortante= ');disp(MySimplexTable(pvt_row+2,1));
    
    end
    BV(pvt_row)=pvt_col;

    
    pivot=A(pvt_row,pvt_col);
    disp("pivot = ");
    disp (pivot);
    
    A(pvt_row,:) = A(pvt_row,:)./pivot;
    for i=1:size(A,1)
        if i~=pvt_row
            A(i,:)=A(i,:)-A(i,pvt_col).*A(pvt_row,:);
        end
    end  
    ZjCj=ZjCj-ZjCj(pvt_col).*A(pvt_row,:);
    
    ZCj=[ZjCj;A];
    tab=array2table(ZCj);
    BFS=zeros(1,size(A,2));
    BFS(end)=sum(BFS.*Cost);
    Current_BFS=array2table(BFS);
    [rows,cols]=size(tab);
    MySimplexTable(pvt_row+2,1)=MySimplexTable(1,pvt_col+1);
    x=[-ZCj(1,size(ZCj,2))];
    y=array2table(x);
    
    for i=1:size(ZCj,1)
        for j=1:size(ZCj,2)
            if (i==1 && j==size(ZCj,2))
                MySimplexTable(i+1,j+1)=y(1,1);
            else
            MySimplexTable(i+1,j+1)=tab(i,j);
            end
        end
    end
        
    disp(MySimplexTable);

else
    RUN=false;
    disp('Solution optimale a ete atteinte');
    
end
end
tab=table2array(MySimplexTable);
for i=1:Ineq+1
    disp(strcat(tab(1+i,1),"=",tab(1+i,end)));
end
else
NbreVariables=input('Donner le nombre des variables');
disp('donner les coefficients du Min dans l ordre');
for i=1:NbreVariables
    b(i,1)=input(strcat('B [', num2str(i),']=   '));
end
Ineq=input('Donner le nombre des inequations');
disp ('Donner les coefficients des varialbles dans les inequations dans l ordre');
for i=1:Ineq
    for j=1:NbreVariables
            Infos(i,j)=input(strcat('Info [', num2str(i),',',num2str(j),']=   '));
    end
end
Info=Infos.';
disp('donner les coefficients à droite de signe des inequations dans l ordre');
for i=1:Ineq
   C(i)=input(strcat('Coeff [', num2str(i),']=   '));
end    
 s=eye(size(Info,1));

A=[Info s b];
K=["Val";"Zj-Cj"];
for i=1:NbreVariables
    K(i+2,1)=strcat("e_",num2str(i));
end
for i=1:Ineq
    J(i)=strcat("y_",num2str(i));
end
for i=Ineq+1:Ineq+NbreVariables
    J(i)=strcat("e_",num2str(i-NbreVariables-1));
end
J(NbreVariables+Ineq+1)="bj";

K1=array2table(K);

Cost=zeros(1,size(A,2));
Cost(1:Ineq)=C;

BV= Ineq+1:1:size(A,2)-1;

ZjCj=Cost-Cost(BV)*A;
ZCj=[ZjCj ; A];
B=[J ; ZCj];
MySimplexTable=array2table(B);
MySimplexTable=[K1 MySimplexTable];
disp(MySimplexTable);
while RUN
if any(ZjCj>0);
     ZC=ZjCj(1:end-1);
       [EnterCol,pvt_col]=max(ZC);

        
     fprintf('Variable entrante est '); disp(MySimplexTable(1,pvt_col+1));
   sol= A(:,end);
    Column=A(:,pvt_col);
    if all(A(:,pvt_col)<=0)
        error('Impossible de resoudre, toutes les valeurs entrees <=0 dans la colonne %d',pvt_col);
    else
    for i=1:size(Column,1)
        if Column(i)>0
            ratio(i)=sol(i)./Column(i);
        else
            ratio(i)=inf;
        end
    end
    end
    [MinRatio,pvt_row]=min(ratio);
    disp('Variable sortante= ');disp(MySimplexTable(pvt_row+2,1));
    
    
    BV(pvt_row)=pvt_col;

    
    pivot=A(pvt_row,pvt_col);
    disp("pivot = ");
    disp (pivot);
    
    A(pvt_row,:) = A(pvt_row,:)./pivot;
    for i=1:size(A,1)
        if i~=pvt_row
            A(i,:)=A(i,:)-A(i,pvt_col).*A(pvt_row,:);
        end
    end  
    ZjCj=ZjCj-ZjCj(pvt_col).*A(pvt_row,:);
    
    ZCj=[ZjCj;A];
    tab=array2table(ZCj);
    BFS=zeros(1,size(A,2));
    BFS(end)=sum(BFS.*Cost);
    Current_BFS=array2table(BFS);
    [rows,cols]=size(tab);
    MySimplexTable(pvt_row+2,1)=MySimplexTable(1,pvt_col+1);
    x=[-ZCj(1,size(ZCj,2))];
    y=array2table(x);
    
    for i=1:size(ZCj,1)
        for j=1:size(ZCj,2)
            if (i==1 && j==size(ZCj,2))
                MySimplexTable(i+1,j+1)=y(1,1);
            else
            MySimplexTable(i+1,j+1)=tab(i,j);
            end
        end
    end
        
    disp(MySimplexTable);

else
    RUN=false;
    disp('Solution optimale a ete atteinte');
end
end
tab=table2array(MySimplexTable);
for i=1:NbreVariables
    disp(strcat("x_",num2str(i),"= "));
    x=-str2num(tab(2,Ineq+1+i));
    disp(x);
end
disp("Min =");
disp(MySimplexTable(2,end));
    
end

 
