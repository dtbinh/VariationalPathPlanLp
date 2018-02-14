function och = VariationalPathPlanFunction()

och.icf = @VariationalPathPlanFunction_icf; 
och.tcf = @VariationalPathPlanFunction_tcf; 
och.fcf = @VariationalPathPlanFunction_fcf; 

och.nlicf = @VariationalPathPlanFunction_nlicf; 
och.nltcf = @VariationalPathPlanFunction_nltcf; 
och.nlfcf = @VariationalPathPlanFunction_nlfcf; 
och.nlgcf = @VariationalPathPlanFunction_nlgcf; 



function [f,df] = VariationalPathPlanFunction_icf(x,xd,xdd,xddd,xdddd,y,yd,ydd,yddd,ydddd,th,thd,thdd,thddd,thdddd,mu,mud,tf)
global qf qdf sigma m J tau obs obsOrder obsScale B1 B2 


end


function [f,df] = VariationalPathPlanFunction_tcf(x,xd,xdd,xddd,xdddd,y,yd,ydd,yddd,ydddd,th,thd,thdd,thddd,thdddd,mu,mud,tf)
global qf qdf sigma m J tau obs obsOrder obsScale B1 B2 


end


function [f,df] = VariationalPathPlanFunction_fcf(x,xd,xdd,xddd,xdddd,y,yd,ydd,yddd,ydddd,th,thd,thdd,thddd,thdddd,mu,mud,tf)
global qf qdf sigma m J tau obs obsOrder obsScale B1 B2 
	 f(1,:) = tf;

	 df(1,:)= zeros(size(x)); 
	 df(2,:)= zeros(size(xd)); 
	 df(3,:)= zeros(size(xdd)); 
	 df(4,:)= zeros(size(xddd)); 
	 df(5,:)= zeros(size(xdddd)); 
	 df(6,:)= zeros(size(y)); 
	 df(7,:)= zeros(size(yd)); 
	 df(8,:)= zeros(size(ydd)); 
	 df(9,:)= zeros(size(yddd)); 
	 df(10,:)= zeros(size(ydddd)); 
	 df(11,:)= zeros(size(th)); 
	 df(12,:)= zeros(size(thd)); 
	 df(13,:)= zeros(size(thdd)); 
	 df(14,:)= zeros(size(thddd)); 
	 df(15,:)= zeros(size(thdddd)); 
	 df(16,:)= zeros(size(mu)); 
	 df(17,:)= zeros(size(mud)); 
	 df(18,:)= 1.*ones(size(tf)); 

end


function [f,df] = VariationalPathPlanFunction_nlicf(x,xd,xdd,xddd,xdddd,y,yd,ydd,yddd,ydddd,th,thd,thdd,thddd,thdddd,mu,mud,tf)
global qf qdf sigma m J tau obs obsOrder obsScale B1 B2 


end


function [f,df] = VariationalPathPlanFunction_nltcf(x,xd,xdd,xddd,xdddd,y,yd,ydd,yddd,ydddd,th,thd,thdd,thddd,thdddd,mu,mud,tf)
global qf qdf sigma m J tau obs obsOrder obsScale B1 B2 
	 f(1,:) = xd.*sin(th)-yd.*cos(th);
	 f(2,:) = (((x-obs(1))./obsScale(1)).^obsOrder+((y-obs(2))./obsScale(2)).^obsOrder).^(1./obsOrder)-1;
	 f(3,:) = B1.^2.*tf.^2-xd.^2-yd.^2;
	 f(4,:) = B2.^2.*tf.^2-thd.^2;
	 f(5,:) = -thdddd+sigma.*thdd.*tf.^2-mu./J.*(xd.*tf.^3.*cos(th)+yd.*tf.^3.*sin(th));
	 f(6,:) = -xdddd+sigma.*xdd.*tf.^2+tau.*obsOrder./2.*((x-obs(1))./obsScale(1)).^(obsOrder-1)./obsScale(1)./(((x-obs(1))./obsScale(1)).^obsOrder+((y-obs(2))./obsScale(2)).^obsOrder-1.^obsOrder).^2+1./m.*(mud.*tf.^3.*sin(th)+mu.*thd.*tf.^3.*cos(th));
	 f(7,:) = -ydddd+sigma.*ydd.*tf.^2+tau.*obsOrder./2.*((y-obs(2))./obsScale(2)).^(obsOrder-1)./obsScale(2)./(((x-obs(1))./obsScale(1)).^obsOrder+((y-obs(2))./obsScale(2)).^obsOrder-1.^obsOrder).^2+1./m.*(-mud.*tf.^3.*cos(th)+mu.*thd.*tf.^3.*sin(th));

	 df(1,1,:) = zeros(size(x)); 
	 df(2,1,:) = sin(th); 
	 df(3,1,:) = zeros(size(xdd)); 
	 df(4,1,:) = zeros(size(xddd)); 
	 df(5,1,:) = zeros(size(xdddd)); 
	 df(6,1,:) = zeros(size(y)); 
	 df(7,1,:) = -cos(th); 
	 df(8,1,:) = zeros(size(ydd)); 
	 df(9,1,:) = zeros(size(yddd)); 
	 df(10,1,:) = zeros(size(ydddd)); 
	 df(11,1,:) = xd.*cos(th)+yd.*sin(th); 
	 df(12,1,:) = zeros(size(thd)); 
	 df(13,1,:) = zeros(size(thdd)); 
	 df(14,1,:) = zeros(size(thddd)); 
	 df(15,1,:) = zeros(size(thdddd)); 
	 df(16,1,:) = zeros(size(mu)); 
	 df(17,1,:) = zeros(size(mud)); 
	 df(18,1,:) = zeros(size(tf)); 
	 df(1,2,:) = (((x-obs(1))./obsScale(1)).^obsOrder+((y-obs(2))./obsScale(2)).^obsOrder).^(1./obsOrder).*((x-obs(1))./obsScale(1)).^obsOrder./(x-obs(1))./(((x-obs(1))./obsScale(1)).^obsOrder+((y-obs(2))./obsScale(2)).^obsOrder); 
	 df(2,2,:) = zeros(size(xd)); 
	 df(3,2,:) = zeros(size(xdd)); 
	 df(4,2,:) = zeros(size(xddd)); 
	 df(5,2,:) = zeros(size(xdddd)); 
	 df(6,2,:) = (((x-obs(1))./obsScale(1)).^obsOrder+((y-obs(2))./obsScale(2)).^obsOrder).^(1./obsOrder).*((y-obs(2))./obsScale(2)).^obsOrder./(y-obs(2))./(((x-obs(1))./obsScale(1)).^obsOrder+((y-obs(2))./obsScale(2)).^obsOrder); 
	 df(7,2,:) = zeros(size(yd)); 
	 df(8,2,:) = zeros(size(ydd)); 
	 df(9,2,:) = zeros(size(yddd)); 
	 df(10,2,:) = zeros(size(ydddd)); 
	 df(11,2,:) = zeros(size(th)); 
	 df(12,2,:) = zeros(size(thd)); 
	 df(13,2,:) = zeros(size(thdd)); 
	 df(14,2,:) = zeros(size(thddd)); 
	 df(15,2,:) = zeros(size(thdddd)); 
	 df(16,2,:) = zeros(size(mu)); 
	 df(17,2,:) = zeros(size(mud)); 
	 df(18,2,:) = zeros(size(tf)); 
	 df(1,3,:) = zeros(size(x)); 
	 df(2,3,:) = -2.*xd; 
	 df(3,3,:) = zeros(size(xdd)); 
	 df(4,3,:) = zeros(size(xddd)); 
	 df(5,3,:) = zeros(size(xdddd)); 
	 df(6,3,:) = zeros(size(y)); 
	 df(7,3,:) = -2.*yd; 
	 df(8,3,:) = zeros(size(ydd)); 
	 df(9,3,:) = zeros(size(yddd)); 
	 df(10,3,:) = zeros(size(ydddd)); 
	 df(11,3,:) = zeros(size(th)); 
	 df(12,3,:) = zeros(size(thd)); 
	 df(13,3,:) = zeros(size(thdd)); 
	 df(14,3,:) = zeros(size(thddd)); 
	 df(15,3,:) = zeros(size(thdddd)); 
	 df(16,3,:) = zeros(size(mu)); 
	 df(17,3,:) = zeros(size(mud)); 
	 df(18,3,:) = 2.*B1.^2.*tf; 
	 df(1,4,:) = zeros(size(x)); 
	 df(2,4,:) = zeros(size(xd)); 
	 df(3,4,:) = zeros(size(xdd)); 
	 df(4,4,:) = zeros(size(xddd)); 
	 df(5,4,:) = zeros(size(xdddd)); 
	 df(6,4,:) = zeros(size(y)); 
	 df(7,4,:) = zeros(size(yd)); 
	 df(8,4,:) = zeros(size(ydd)); 
	 df(9,4,:) = zeros(size(yddd)); 
	 df(10,4,:) = zeros(size(ydddd)); 
	 df(11,4,:) = zeros(size(th)); 
	 df(12,4,:) = -2.*thd; 
	 df(13,4,:) = zeros(size(thdd)); 
	 df(14,4,:) = zeros(size(thddd)); 
	 df(15,4,:) = zeros(size(thdddd)); 
	 df(16,4,:) = zeros(size(mu)); 
	 df(17,4,:) = zeros(size(mud)); 
	 df(18,4,:) = 2.*B2.^2.*tf; 
	 df(1,5,:) = zeros(size(x)); 
	 df(2,5,:) = -mu./J.*tf.^3.*cos(th); 
	 df(3,5,:) = zeros(size(xdd)); 
	 df(4,5,:) = zeros(size(xddd)); 
	 df(5,5,:) = zeros(size(xdddd)); 
	 df(6,5,:) = zeros(size(y)); 
	 df(7,5,:) = -mu./J.*tf.^3.*sin(th); 
	 df(8,5,:) = zeros(size(ydd)); 
	 df(9,5,:) = zeros(size(yddd)); 
	 df(10,5,:) = zeros(size(ydddd)); 
	 df(11,5,:) = -mu./J.*(-xd.*tf.^3.*sin(th)+yd.*tf.^3.*cos(th)); 
	 df(12,5,:) = zeros(size(thd)); 
	 df(13,5,:) = sigma.*tf.^2; 
	 df(14,5,:) = zeros(size(thddd)); 
	 df(15,5,:) = -1.*ones(size(thdddd)); 
	 df(16,5,:) = -1./J.*(xd.*tf.^3.*cos(th)+yd.*tf.^3.*sin(th)); 
	 df(17,5,:) = zeros(size(mud)); 
	 df(18,5,:) = 2.*sigma.*thdd.*tf-mu./J.*(3.*xd.*tf.^2.*cos(th)+3.*yd.*tf.^2.*sin(th)); 
	 df(1,6,:) = 1./2.*tau.*obsOrder.*((x-obs(1))./obsScale(1)).^(obsOrder-1).*(obsOrder-1)./(x-obs(1))./obsScale(1)./(((x-obs(1))./obsScale(1)).^obsOrder+((y-obs(2))./obsScale(2)).^obsOrder-1).^2-tau.*obsOrder.^2.*((x-obs(1))./obsScale(1)).^(obsOrder-1)./obsScale(1)./(((x-obs(1))./obsScale(1)).^obsOrder+((y-obs(2))./obsScale(2)).^obsOrder-1).^3.*((x-obs(1))./obsScale(1)).^obsOrder./(x-obs(1)); 
	 df(2,6,:) = zeros(size(xd)); 
	 df(3,6,:) = sigma.*tf.^2; 
	 df(4,6,:) = zeros(size(xddd)); 
	 df(5,6,:) = -1.*ones(size(xdddd)); 
	 df(6,6,:) = -tau.*obsOrder.^2.*((x-obs(1))./obsScale(1)).^(obsOrder-1)./obsScale(1)./(((x-obs(1))./obsScale(1)).^obsOrder+((y-obs(2))./obsScale(2)).^obsOrder-1).^3.*((y-obs(2))./obsScale(2)).^obsOrder./(y-obs(2)); 
	 df(7,6,:) = zeros(size(yd)); 
	 df(8,6,:) = zeros(size(ydd)); 
	 df(9,6,:) = zeros(size(yddd)); 
	 df(10,6,:) = zeros(size(ydddd)); 
	 df(11,6,:) = 1./m.*(mud.*tf.^3.*cos(th)-mu.*thd.*tf.^3.*sin(th)); 
	 df(12,6,:) = 1./m.*mu.*tf.^3.*cos(th); 
	 df(13,6,:) = zeros(size(thdd)); 
	 df(14,6,:) = zeros(size(thddd)); 
	 df(15,6,:) = zeros(size(thdddd)); 
	 df(16,6,:) = 1./m.*thd.*tf.^3.*cos(th); 
	 df(17,6,:) = 1./m.*tf.^3.*sin(th); 
	 df(18,6,:) = 2.*sigma.*xdd.*tf+1./m.*(3.*mud.*tf.^2.*sin(th)+3.*mu.*thd.*tf.^2.*cos(th)); 
	 df(1,7,:) = -tau.*obsOrder.^2.*((y-obs(2))./obsScale(2)).^(obsOrder-1)./obsScale(2)./(((x-obs(1))./obsScale(1)).^obsOrder+((y-obs(2))./obsScale(2)).^obsOrder-1).^3.*((x-obs(1))./obsScale(1)).^obsOrder./(x-obs(1)); 
	 df(2,7,:) = zeros(size(xd)); 
	 df(3,7,:) = zeros(size(xdd)); 
	 df(4,7,:) = zeros(size(xddd)); 
	 df(5,7,:) = zeros(size(xdddd)); 
	 df(6,7,:) = 1./2.*tau.*obsOrder.*((y-obs(2))./obsScale(2)).^(obsOrder-1).*(obsOrder-1)./(y-obs(2))./obsScale(2)./(((x-obs(1))./obsScale(1)).^obsOrder+((y-obs(2))./obsScale(2)).^obsOrder-1).^2-tau.*obsOrder.^2.*((y-obs(2))./obsScale(2)).^(obsOrder-1)./obsScale(2)./(((x-obs(1))./obsScale(1)).^obsOrder+((y-obs(2))./obsScale(2)).^obsOrder-1).^3.*((y-obs(2))./obsScale(2)).^obsOrder./(y-obs(2)); 
	 df(7,7,:) = zeros(size(yd)); 
	 df(8,7,:) = sigma.*tf.^2; 
	 df(9,7,:) = zeros(size(yddd)); 
	 df(10,7,:) = -1.*ones(size(ydddd)); 
	 df(11,7,:) = 1./m.*(mud.*tf.^3.*sin(th)+mu.*thd.*tf.^3.*cos(th)); 
	 df(12,7,:) = 1./m.*mu.*tf.^3.*sin(th); 
	 df(13,7,:) = zeros(size(thdd)); 
	 df(14,7,:) = zeros(size(thddd)); 
	 df(15,7,:) = zeros(size(thdddd)); 
	 df(16,7,:) = 1./m.*thd.*tf.^3.*sin(th); 
	 df(17,7,:) = -1./m.*tf.^3.*cos(th); 
	 df(18,7,:) = 2.*sigma.*ydd.*tf+1./m.*(-3.*mud.*tf.^2.*cos(th)+3.*mu.*thd.*tf.^2.*sin(th)); 

end


function [f,df] = VariationalPathPlanFunction_nlfcf(x,xd,xdd,xddd,xdddd,y,yd,ydd,yddd,ydddd,th,thd,thdd,thddd,thdddd,mu,mud,tf)
global qf qdf sigma m J tau obs obsOrder obsScale B1 B2 


end


function [f,df] = VariationalPathPlanFunction_nlgcf(x,xd,xdd,xddd,xdddd,y,yd,ydd,yddd,ydddd,th,thd,thdd,thddd,thdddd,mu,mud,tf)
global qf qdf sigma m J tau obs obsOrder obsScale B1 B2 


end


end