function och = VariationalPathPlan()

och.icf = @VariationalPathPlan_icf; 
och.tcf = @VariationalPathPlan_tcf; 
och.fcf = @VariationalPathPlan_fcf; 

och.nlicf = @VariationalPathPlan_nlicf; 
och.nltcf = @VariationalPathPlan_nltcf; 
och.nlfcf = @VariationalPathPlan_nlfcf; 
och.nlgcf = @VariationalPathPlan_nlgcf; 



function [f,df] = VariationalPathPlan_icf(x,xd,xdd,,y,yd,ydd,,th,thd,thdd,,mu,mud)
global qf qdf sigma m J tau 


end


function [f,df] = VariationalPathPlan_tcf(x,xd,xdd,,y,yd,ydd,,th,thd,thdd,,mu,mud)
global qf qdf sigma m J tau 


end


function [f,df] = VariationalPathPlan_fcf(x,xd,xdd,,y,yd,ydd,,th,thd,thdd,,mu,mud)
global qf qdf sigma m J tau 
	 f(1,:) = (x-qf(1)).^2+(y-qf(2)).^2+(th-qf(3)).^2;

	 df(1,:)= 2.*x-2.*qf(1); 
	 df(2,:)= zeros(size(xd)); 
	 df(3,:)= zeros(size(xdd)); 
	 df(4,:)= zeros(size()); 
	 df(5,:)= 2.*y-2.*qf(2); 
	 df(6,:)= zeros(size(yd)); 
	 df(7,:)= zeros(size(ydd)); 
	 df(8,:)= zeros(size()); 
	 df(9,:)= 2.*th-2.*qf(3); 
	 df(10,:)= zeros(size(thd)); 
	 df(11,:)= zeros(size(thdd)); 
	 df(12,:)= zeros(size()); 
	 df(13,:)= zeros(size(mu)); 
	 df(14,:)= zeros(size(mud)); 

end


function [f,df] = VariationalPathPlan_nlicf(x,xd,xdd,,y,yd,ydd,,th,thd,thdd,,mu,mud)
global qf qdf sigma m J tau 


end


function [f,df] = VariationalPathPlan_nltcf(x,xd,xdd,,y,yd,ydd,,th,thd,thdd,,mu,mud)
global qf qdf sigma m J tau 
	 f(1,:) = xd.*sin(th)-yd.*cos(th);
	 f(2,:) = -thdddd+sigma.*thdd-mu./J.*(xd.*cos(th)+yd.*sin(th));
	 f(3,:) = -xdddd+sigma.*xdd+tau.*x./((x-obs(1)).^2+(y-obs(2)).^2-obs(3).^2)+1./m.*(mud.*sin(th)+mu.*thd.*cos(th));
	 f(4,:) = -ydddd+sigma.*ydd+tau.*y./((x-obs(1)).^2+(y-obs(2)).^2-obs(3).^2)+1./m.*(-mud.*cos(th)+mu.*thd.*sin(th));

	 df(1,1,:) = zeros(size(x)); 
	 df(2,1,:) = sin(th); 
	 df(3,1,:) = zeros(size(xdd)); 
	 df(4,1,:) = zeros(size()); 
	 df(5,1,:) = zeros(size(y)); 
	 df(6,1,:) = -cos(th); 
	 df(7,1,:) = zeros(size(ydd)); 
	 df(8,1,:) = zeros(size()); 
	 df(9,1,:) = xd.*cos(th)+yd.*sin(th); 
	 df(10,1,:) = zeros(size(thd)); 
	 df(11,1,:) = zeros(size(thdd)); 
	 df(12,1,:) = zeros(size()); 
	 df(13,1,:) = zeros(size(mu)); 
	 df(14,1,:) = zeros(size(mud)); 
	 df(1,2,:) = zeros(size(x)); 
	 df(2,2,:) = -mu./J.*cos(th); 
	 df(3,2,:) = zeros(size(xdd)); 
	 df(4,2,:) = zeros(size()); 
	 df(5,2,:) = zeros(size(y)); 
	 df(6,2,:) = -mu./J.*sin(th); 
	 df(7,2,:) = zeros(size(ydd)); 
	 df(8,2,:) = zeros(size()); 
	 df(9,2,:) = -mu./J.*(-xd.*sin(th)+yd.*cos(th)); 
	 df(10,2,:) = zeros(size(thd)); 
	 df(11,2,:) = sigma; 
	 df(12,2,:) = zeros(size()); 
	 df(13,2,:) = -1./J.*(xd.*cos(th)+yd.*sin(th)); 
	 df(14,2,:) = zeros(size(mud)); 
	 df(1,3,:) = tau./((x-obs(1)).^2+(y-obs(2)).^2-obs(3).^2)-tau.*x./((x-obs(1)).^2+(y-obs(2)).^2-obs(3).^2).^2.*(2.*x-2.*obs(1)); 
	 df(2,3,:) = zeros(size(xd)); 
	 df(3,3,:) = sigma; 
	 df(4,3,:) = zeros(size()); 
	 df(5,3,:) = -tau.*x./((x-obs(1)).^2+(y-obs(2)).^2-obs(3).^2).^2.*(2.*y-2.*obs(2)); 
	 df(6,3,:) = zeros(size(yd)); 
	 df(7,3,:) = zeros(size(ydd)); 
	 df(8,3,:) = zeros(size()); 
	 df(9,3,:) = 1./m.*(mud.*cos(th)-mu.*thd.*sin(th)); 
	 df(10,3,:) = 1./m.*mu.*cos(th); 
	 df(11,3,:) = zeros(size(thdd)); 
	 df(12,3,:) = zeros(size()); 
	 df(13,3,:) = 1./m.*thd.*cos(th); 
	 df(14,3,:) = 1./m.*sin(th); 
	 df(1,4,:) = -tau.*y./((x-obs(1)).^2+(y-obs(2)).^2-obs(3).^2).^2.*(2.*x-2.*obs(1)); 
	 df(2,4,:) = zeros(size(xd)); 
	 df(3,4,:) = zeros(size(xdd)); 
	 df(4,4,:) = zeros(size()); 
	 df(5,4,:) = tau./((x-obs(1)).^2+(y-obs(2)).^2-obs(3).^2)-tau.*y./((x-obs(1)).^2+(y-obs(2)).^2-obs(3).^2).^2.*(2.*y-2.*obs(2)); 
	 df(6,4,:) = zeros(size(yd)); 
	 df(7,4,:) = sigma; 
	 df(8,4,:) = zeros(size()); 
	 df(9,4,:) = 1./m.*(mud.*sin(th)+mu.*thd.*cos(th)); 
	 df(10,4,:) = 1./m.*mu.*sin(th); 
	 df(11,4,:) = zeros(size(thdd)); 
	 df(12,4,:) = zeros(size()); 
	 df(13,4,:) = 1./m.*thd.*sin(th); 
	 df(14,4,:) = -1./m.*cos(th); 

end


function [f,df] = VariationalPathPlan_nlfcf(x,xd,xdd,,y,yd,ydd,,th,thd,thdd,,mu,mud)
global qf qdf sigma m J tau 


end


function [f,df] = VariationalPathPlan_nlgcf(x,xd,xdd,,y,yd,ydd,,th,thd,thdd,,mu,mud)
global qf qdf sigma m J tau 


end


end