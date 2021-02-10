function [t_measstop,t_measstart] = parsemeas(t,meas)
% mat = [t meas];
% meas = [];
% for i = 1:numel(t)
%     if (~isnan(mat(i,2))||~isnan(mat(i,3))||~isnan(mat(i,4)))
%         meas = [meas;mat(i,:)];
%     end
% end
% out = meas;
% counter = 0;
% for i = 1:size(meas,1)
%     for j = 2:4
%         if isnan(meas(i,j)) ~= 1
%             counter = counter + 1;
%         end
%     end
% end
t_measstop = [];
t_measstart = [];
for i = 2:numel(t)
    if isnan(meas(i,:))
        if ~all(isnan(meas(i-1,:)))
        t_measstop = [t_measstop i-1];
        end
    end
    
    if ~all(isnan(meas(i,:)))
        if isnan(meas(i-1,:))
        t_measstart = [t_measstart i];
        end
    end
end