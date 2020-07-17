function [X,data_2] = argmin_X_full(Z,U,rho,data,A,II,get_energy,usemex,...
    aggregate,opnorm,weights,NN,BXX,BXY,BYX,BYY,BX,BY,C)

X = zeros(size(A,1),1);
AZ = A*Z;

if aggregate==1
    if opnorm
        for i=0:(length(II)-1)
            AZi = [AZ(4*i+1),AZ(4*i+2);AZ(4*i+3),AZ(4*i+4)];
            MM = AZi-[U(4*i+1),U(4*i+2);U(4*i+3),U(4*i+4)];
            [T,S,V] = svd(MM);
            h1 = S(1,1);
            if (S(2,2)-(1/rho))>0
                h2 = S(2,2)-(1/rho);
            else
                h2 = 0;
            end
            Hval = T*[h1,0;0,h2]*(V');
            X(4*i+1) = Hval(1,1);
            X(4*i+2) = Hval(1,2);
            X(4*i+3) = Hval(2,1);
            X(4*i+4) = Hval(2,2);
            
        end
    else
        
        if usemex
            
            % weights = ones(length(II),1);
            X = argmin_X_full_mex(AZ,U,II,rho,weights);
            
        else
            
            
            
            for i=0:((size(A,1)/4)-1)
                % disp(i)
                %     yalmip("clear")
                %     H = sdpvar(2,2,'full');
                %     AZi = [AZ(4*i+1);AZ(4*i+2);AZ(4*i+3);AZ(4*i+4)];
                %     objective = norm(H,'nuclear')+(rho/2)*norm([H(1,1);H(1,2);H(2,1);H(2,2)] - ...
                %         AZi+[U(4*i+1);U(4*i+2);U(4*i+3);U(4*i+4)])^2;
                %     options = sdpsettings('verbose',0,'debug',1,'solver','mosek',...
                %         'cachesolvers',1);
                %     constraints = [];
                %     sol = optimize(constraints,objective,options);
                %     Hval = value(H);
                %     X(4*i+1) = Hval(1,1);
                %     X(4*i+2) = Hval(1,2);
                %     X(4*i+3) = Hval(2,1);
                %     X(4*i+4) = Hval(2,2);
                
                AZi = [AZ(4*i+1),AZ(4*i+2);AZ(4*i+3),AZ(4*i+4)];
                MM = AZi-[U(4*i+1),U(4*i+2);U(4*i+3),U(4*i+4)];
                [T,S,V] = svd(MM);
                
                % if S(1,1) == 0
                %     h1 = 0;
                % else
                %     h1 = S(1,1)-(sign(S(1,1))/rho);
                % end
                rho_with_weight = rho/weights(i+1);
                
                
                if S(1,1)-(1/rho_with_weight)>0
                    h1 = S(1,1)-(1/rho_with_weight);
                elseif S(1,1)+(1/rho_with_weight)<0
                    h1 = S(1,1)+(1/rho_with_weight);
                else
                    h1 = 0;
                end
                
                if S(2,2)-(1/rho_with_weight)>0
                    h2 = S(2,2)-(1/rho_with_weight);
                elseif S(2,2)+(1/rho_with_weight)<0
                    h2 = S(2,2)+(1/rho_with_weight);
                else
                    h2 = 0;
                end
                
                Hval = T*[h1,0;0,h2]*(V');
                X(4*i+1) = Hval(1,1);
                X(4*i+2) = Hval(1,2);
                X(4*i+3) = Hval(2,1);
                X(4*i+4) = Hval(2,2);
            end
            
        end
        
    end
elseif aggregate==2
    if opnorm
        for i=0:(length(II)-1)
            AZi = [AZ(4*i+1),AZ(4*i+2);AZ(4*i+3),AZ(4*i+4)];
            MM = AZi-[U(4*i+1),U(4*i+2);U(4*i+3),U(4*i+4)];
            [T,S,V] = svd(MM);
            h1 = S(1,1);
            h2 = (rho/(rho+2))*S(2,2);
            Hval = T*[h1,0;0,h2]*(V');
            X(4*i+1) = Hval(1,1);
            X(4*i+2) = Hval(1,2);
            X(4*i+3) = Hval(2,1);
            X(4*i+4) = Hval(2,2);
            
        end
    else
        
        if usemex
            X = argmin_X_full_L2_mex(AZ,U,II,rho);
        else
            
            for i=0:(length(II)-1)
                
                AZi = [AZ(4*i+1),AZ(4*i+2);AZ(4*i+3),AZ(4*i+4)];
                MM = AZi-[U(4*i+1),U(4*i+2);U(4*i+3),U(4*i+4)];
                [T,S,V] = svd(MM);
                
                if (((rho+2)*S(1,1))-(2*S(2,2)))>0 && (((rho+2)*S(2,2))-(2*S(1,1)))>0
                    h1 = (((rho+2)*S(1,1))-(2*S(2,2)))/(4+rho);
                    h2 = (((rho+2)*S(2,2))-(2*S(1,1)))/(4+rho);
                else
                    assert(S(1,1)>=S(2,2));
                    h1 = (rho/(rho+2))*S(1,1);
                    h2 = 0;
                end
                
                
                
                
                Hval = T*[h1,0;0,h2]*(V');
                X(4*i+1) = Hval(1,1);
                X(4*i+2) = Hval(1,2);
                X(4*i+3) = Hval(2,1);
                X(4*i+4) = Hval(2,2);
                
            end
            
        end
    end
elseif aggregate==3
    for i=0:(length(II)-1)
        
        AZi = [AZ(4*i+1),AZ(4*i+2);AZ(4*i+3),AZ(4*i+4)];
        MM = AZi-[U(4*i+1),U(4*i+2);U(4*i+3),U(4*i+4)];
        [T,S,V] = svd(MM);
        
        
        
        h1 = (rho/(rho+2))*S(1,1);
        h2 = (rho/(rho+2))*S(2,2);
        
        
        
        
        
        Hval = T*[h1,0;0,h2]*(V');
        X(4*i+1) = Hval(1,1);
        X(4*i+2) = Hval(1,2);
        X(4*i+3) = Hval(2,1);
        X(4*i+4) = Hval(2,2);
        
    end
    
end


data_2 = data;
if get_energy
    data_2.energy = 0;
    for i=0:(length(II)-1)
        S = svd([X(4*i+1),X(4*i+2);X(4*i+3),X(4*i+4)]);
        data_2.energy = data_2.energy + sum(S);
    end
end


% for i=1:length(X)
%     % assert(X(i)==Xtest(i))
%     disp(X(i))
%     disp(Xtest(i))
%     pause
% end




end

