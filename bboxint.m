function intarea = bboxint(bbox1,bbox2)
%function intarea = bboxint(bbox1,bbox2)
%
% Return area of intersection of two bounding boxes.
% CALLS: BBOX2RECT, RECTINT
% 
% Last Saved Time-stamp: <Thu 2016-05-19 12:36:43 Eastern Daylight Time lew.gramer>

  intarea = rectint(bbox2rect(bbox1),bbox2rect(bbox2));

return;

