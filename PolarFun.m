classdef PolarFun
    methods (Static)
        
        function Tval = EvaluoT(T,r,alpha)
            Tval = eval(T);
        end
        
        function [rs,alphas] = GeneraRsAlphas(Ncom,rmean,alphamean,ancho,alto)
            
            rdev = .5;
            alphadev = .1;
            %ancho = 0.2;
%             alto = 1.25;            
            alphas = nan(1,Ncom);
            rs = nan(1,Ncom);
            
            xv1 = linspace(rmean-ancho,rmean+ancho,10);
            xv2 = fliplr(xv1);
            limup = rmean*2/alphamean+alto/2;
            limdown = rmean*2/alphamean-alto/2;
            yv = [limdown./xv1 limup./xv2 limdown/xv1(1)];
            xv = [xv1 xv2 xv1(1)];
            for indp = 1:Ncom
                good = 0;
                red = 1;
                while good == 0
                    red = rdev*randn+rmean;
                    alphaed = alphadev*randn+alphamean;
                    good = inpolygon(red,2./alphaed,xv,yv);
                end
                alphas(indp) = alphaed;
                rs(indp) = red;
            end
            
        end
                
        function dydt = Acoplado(t,y,rs,alphas,mytrend,years,mlf1,mlf2,mlf3,C)
            
            Nwords = length(C);
            
            dydt = zeros(Nwords*3,1);
            
            radios = y(1+(0:(Nwords-1))*3);
            titas = y(2+(0:(Nwords-1))*3);
            ws = y(3+(0:(Nwords-1))*3);
            
            for indp = 1:Nwords
                [k,kpunto] = PolarFun.myinterp1(years,mytrend(indp,:),t);
                
                radio = radios(indp);
                tita = titas(indp);
                w = ws(indp);
                
                Ncom = sum(C(indp,:)~=0)+1;
                
                r = rs(indp);
                alpha = alphas(indp);
                
                dydt(1+(indp-1)*3) = mlf1(alpha,k,kpunto,r,radio,tita,w);
                dydt(2+(indp-1)*3) = mlf2(alpha,k,kpunto,r,radio,tita,w)...
                    +sum(C(indp,:)'.*(sin(titas-tita)))/Ncom;
                dydt(3+(indp-1)*3) = mlf3(alpha,k,kpunto,r,radio,tita,w);
            end

        end
        
        function dydt = AcopladoDF(t,y,rs,alphas,mytrend,years,mlf1,mlf2,mlf3,CD,CF,ND,NF,tstart,indac,cuantosac)
            cc = clock;
            fprintf('%i/%i t = %1.5f time = %1.2f s %02g:%02g:%02g\n',indac,...
                cuantosac,t,toc(tstart),cc(4),cc(5),round(cc(6)))
            Nwords = length(ND);
            
            dydt = zeros(Nwords*3,1);
            
            radios = y(1+(0:(Nwords-1))*3);
            titas = y(2+(0:(Nwords-1))*3);
            ws = y(3+(0:(Nwords-1))*3);
            for indp = 1:Nwords
                [k,kpunto] = PolarFun.myinterp1(years,mytrend(indp,:),t);
                
                radio = radios(indp);
                tita = titas(indp);
                w = ws(indp);
                                
                r = rs(indp);
                alpha = alphas(indp);
                
                dydt(1+(indp-1)*3) = mlf1(alpha,k,kpunto,r,radio,tita,w);
                dydt(2+(indp-1)*3) = mlf2(alpha,k,kpunto,r,radio,tita,w)...
                    +sum(CD(indp,:)'.*(sin(titas-tita)))/ND(indp) ...
                    +sum(CF(indp,:)'.*(sin(titas-tita)))/NF(indp);
                dydt(3+(indp-1)*3) = mlf3(alpha,k,kpunto,r,radio,tita,w);
            end
        end

        function dydt = Desacoplado(t,y,mytrend,alpha,r,years)
            warning('off','all')
            dydt = zeros(3,1);
            dydt(1) = r*y(1)*(1-y(3)/(PolarFun.myinterp1(years,mytrend,t)));
            dydt(2) = alpha*(y(1)-y(2));
            dydt(3) = alpha*(y(2)-y(3));
        end
        
        function out = EvaluoMLF(mlf,alpha,k,kpunto,r,radio,tita,w)
            out = mlf(alpha,k,kpunto,r,radio,tita,w);
        end
        
        function [trendy,y0] = SacoTrend(y0)
            
            numpeaks = zeros(1,21);
            x = reshape(y0,length(y0),1);
            for count = 1:length(numpeaks)
                omegatrend = 0.0195-0.0005*(count-1);
                trendx = autotrend(x,floor(length(x)/2),omegatrend);
                pks = findpeaks(trendx);
                numpeaks(count) = length(pks);
            end
            
            [~,mindex] = min(numpeaks);
            omegatrend = 0.0195-0.0005*(mindex-1);
            trendy = autotrend(x,floor(length(x)/2),omegatrend);
            trendy(isnan(trendy)) = 0;
            trendy = reshape(trendy,size(y0));
            
        end
        
        function [par,ang] = ParOrden(osc)
            Ncom = size(osc,1);
            angulo = angle(hilbert(osc'));
            rhox = sum(cos(angulo),2);
            rhoy = sum(sin(angulo),2);
            [ang,rho] = cart2pol(rhox,rhoy);
            par = rho/Ncom;
        end
        
        function x = Pol2Real(radio,tita,w,mytrend,rs,alphas,T)
            a = radio.*cos(tita);
            b = radio.*sin(tita);
            u = a+1i*b;
            v = a-1i*b;
            
            Tval = PolarFun.EvaluoT(T,rs,alphas);
            x = real(Tval(1,:)*[u;v;w]+mytrend);
        end
        
        function [cipols,ci,azar] = CiPol(T,rs,alphas,mytrend)
            azar = [rand();rand();rand()];
            ci = mytrend(1)*(1+0.6*(azar-.5));
            Tval = PolarFun.EvaluoT(T,rs,alphas);
            T1 = inv(Tval);
            ciuvw = T1*(ci-mytrend(1));
            cipols = [abs(ciuvw(1));...
                atan2(imag(ciuvw(1)),real(ciuvw(1)));...
                ciuvw(3)];
        end
        
        function [k,kpunto] = myinterp1(x,y,xq)
            
            if xq<x(1) || xq>x(end)
                k = NaN;
                kpunto = NaN;
                return
            end
            
            epsilon = 0.000001;
            
            if xq <= x(1)+epsilon
                y1 = y(1);
                [~,kpunto] = interp1p(x,y,xq+2*epsilon);
                k = y1;
            elseif xq >= x(end)-epsilon
                [~,kpunto] = interp1p(x,y,xq-2*epsilon);
                y2 = y(end);
                k = y2;
            else
                y1 = interp1p(x,y,xq-epsilon);
                y2 = interp1p(x,y,xq+epsilon);
                k = (y2+y1)/2;
                kpunto = (y2-y1)/(2*epsilon);
            end
            
            function [yq, pend] = interp1p(x,y,xq)
                inddesde = find(x<xq,1,'last');
                indhasta = find(x>xq,1,'first');
                pend = (y(indhasta)-y(inddesde))/(x(indhasta)-x(inddesde));
                yq = y(inddesde) + pend*(xq-x(inddesde));
            end
            
            
        end
        
        function [value,isterminal,position] = MyEventFun(t,y,Ncom,years,mytrends,rs,alphas,T,x)
            
            value = zeros(1,Ncom);
            isterminal = ones(1,Ncom);
            position = zeros(1,Ncom);
            radios = y(1+(0:(Ncom-1))*3);
            titas = y(2+(0:(Ncom-1))*3);
            ws = y(3+(0:(Ncom-1))*3);
            for indp = 1:Ncom
                r = rs(indp);
                alpha = alphas(indp);
                k = interp1(years,mytrends(indp,:),t);
                
                radio = radios(indp);
                tita = titas(indp);
                w = ws(indp);
                
                value(indp) = PolarFun.Xeval(x,radio,tita,w,T,k,r,alpha);
            end
        end
        
        function xval = Xeval(x,radio,tita,w,T,k,r,alpha)
            T = eval(T);
            xval = eval(x);
        end
        
        function C = ArmoConexionFD(Ncom,debil,fuerte,i_pal,i_coinciden,tipo)
            Conex = ones(Ncom)-eye(Ncom);
            Conex = Conex*debil;
            tipo = lower(tipo);
            
            if strcmp(tipo,'exp')
                for i = 1:length(i_coinciden)
                    i_fuerte = i_coinciden{i};
                    %                     i_fuerte = i_fuerte{1};
                    i_fuerte_lock = find(ismember(i_pal,i_fuerte));
                    Conex(i_fuerte_lock,i_fuerte_lock) = fuerte;
                    Conex = Conex-diag(diag(Conex));
                end
                C = Conex;
            elseif strcmp(tipo,'rand')
                i_fuerte_lock = i_coinciden;
%                     %                     i_fuerte = i_fuerte{1};
%                     i_fuerte_lock = find(ismember(i_pal,i_fuerte));
                Conex(i_fuerte_lock,i_fuerte_lock) = fuerte;
                Conex = Conex-diag(diag(Conex));
                C = Conex;
            elseif strcmp(tipo,'subcom')
                for i = 1:length(i_coinciden)
                    i_fuerte_lock = i_coinciden{i};
                    %                     i_fuerte = i_fuerte{1};
%                     i_fuerte_lock = find(ismember(i_pal,i_fuerte));
                    Conex(i_fuerte_lock,i_fuerte_lock) = fuerte;
                    Conex = Conex-diag(diag(Conex));
                end
                C = Conex;
            end
            
        end
        
    end
end