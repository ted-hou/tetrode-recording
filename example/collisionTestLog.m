% Mouse, date, channel (TR label), unit, collision, clear response, any (including multiunit) response, notes
collisionTestUnits = {...
    'desmond15', '20200225', 1, 3, false, false, true, 'Complex stim response, no collisions, low SNR unit';...
    'desmond15', '20200225', 2, 1, false, true, true, 'Stim response, no collisions, low SNR unit';...
    'desmond15', '20200225', 2, 3, false, true, true, 'Stim response, no collisions, high SNR unit';...
    'desmond15', '20200225', 7, 1, false, false, true, 'Complex stim response, no collisions, low SNR';...
    'desmond15', '20200226', 1, 1, false, true, true, 'Rare but low latency (1ms) stim response, no collisions, low SNR';...
    'desmond15', '20200226', 2, 1, false, true, true, 'Stim response, no spontaneous response, high SNR unit. !!!CT mostly misaligned!!!!!!!!!!!!!!!!!!';...
    'desmond15', '20200226', 3, 1, false, false, false, 'Ignore. Low SNR.';...
    'desmond15', '20200226', 13, 1, false, false, false, 'High SNR. No stim response.';...
    'desmond15', '20200226', 14, 1, false, false, false, 'High SNR. No stim response.';...
    'desmond15', '20200226', 16, 1, false, false, false, 'Low SNR. No stim response.';...
    'desmond15', '20200227', 1, 1, false, true, true, 'Stim response no collisions.';...
    'desmond15', '20200227', 2, 1, false, false, false, 'Low SNR, no response.';...
    'desmond15', '20200227', 4, 1, false, false, false, 'Low SNR, no response.';...
    'desmond15', '20200227', 9, 1, false, false, false, 'High SNR, no response.';...
    'desmond15', '20200227', 11, 1, true, true, true, 'High SNR, Collisions. Every. Time.';...
    'desmond15', '20200227', 12, 1, true, true, true, 'High SNR, Collisions. Every. Time.';...
    'desmond15', '20200302_2', 1, 1, false, false, false, 'No response.';...
    'desmond15', '20200302_2', 4, 1, false, false, false, 'No response.';...
    'desmond15', '20200302_2', 5, 1, true, true, true, 'High SNR, Collisions. Every. Time.';...
    'desmond15', '20200303', 1, 1, false, false, true, 'Low SNR. Weird response.';...
    'desmond15', '20200303', 2, 1, false, false, false, 'High SNR. No response.';...
    'desmond15', '20200303', 5, 1, false, false, false, 'Low SNR. No response.';...
    'desmond16', '20200225', 7, 1, false, false, false, 'Med SNR, No response.';...
    'desmond16', '20200225', 9, 1, false, false, false, 'Low SNR, No response.';...
    'desmond16', '20200227', 3, 1, false, false, false, 'No response.';...
    'desmond16', '20200227', 5, 1, false, false, false, 'No response.';...
    'desmond16', '20200227', 8, 1, false, false, false, 'No response.';...
    'desmond16', '20200227', 8, 1, false, false, false, 'No response.';...
    'desmond16', '20200228', 2, 1, false, false, false, 'High SNR No response.';...
    'desmond16', '20200228', 5, 1, false, false, true, 'Low SNR, response, no collision.';...
    'desmond16', '20200304', 3, 1, false, false, false, 'High SNR, no response, no collision.';...
    'desmond16', '20200305', 1, 1, false, false, true, 'Low SNR, complex response.';...
    'desmond16', '20200305', 4, 1, false, false, false, 'High SNR, no response.';...
    'desmond16', '20200305', 5, 1, false, false, false, 'High SNR, no response.';...
    'desmond16', '20200305', 6, 1, false, false, false, 'High SNR, no response.';...
    'desmond16', '20200305', 7, 1, false, false, true, 'Med SNR, complex response.';...
    'desmond16', '20200306', 1, 1, false, false, false, 'No response.';...
    'desmond16', '20200306', 3, 1, false, false, false, 'No response.';...
    'desmond16', '20200309', 1, 1, false, false, false, 'No response.';...
    'desmond16', '20200309', 2, 1, false, false, false, 'No response.';...
    'desmond16', '20200309', 6, 1, false, false, false, 'No response.';...
    'desmond16', '20200309', 8, 1, false, false, false, 'No response.';...
    'desmond16', '20200309', 10, 1, false, false, true, 'Complex response.';...
    'desmond16', '20200310', 9, 1, false, false, true, 'Response, no collision';...
    'desmond16', '20200310', 8, 1, false, false, true, 'Complex response, no collision';...
    'desmond16', '20200310', 7, 1, false, false, true, 'Complex response, no collision';...
    'desmond16', '20200310', 6, 1, false, false, true, 'Complex response but no unit response, huge SNR, no collision';...
    'desmond16', '20200310', 5, 1, false, false, false, 'No unit response, no collision';...
    'desmond16', '20200310', 4, 1, false, false, false, 'No unit response, no collision';...
    'desmond16', '20200310', 3, 1, false, false, false, 'No unit response, no collision';...
    'desmond16', '20200310', 2, 1, false, false, false, 'No unit response, no collision';...
    'desmond16', '20200310', 1, 1, false, false, false, 'No unit response, no collision';...
    'desmond16', '20200311', 8, 1, false, false, true, 'Complex response no collision high SNR';...
    'desmond16', '20200311', 7, 1, false, false, true, 'Complex response no collision high SNR';...
    'desmond16', '20200311', 6, 1, false, false, false, 'No response huge SNR';...
    'desmond16', '20200311', 5, 1, false, false, false, 'No response no collision huge SNR';...
    'desmond16', '20200311', 4, 1, false, false, false, 'No response no collision huge SNR';...
    'desmond16', '20200311', 1, 1, false, false, false, 'No response no collision big SNR';...
    'desmond16', '20200312', 15, 1, false, false, false, 'No response';...
    'desmond16', '20200312', 13, 1, false, false, false, 'No response big SNR';...
    'desmond16', '20200312', 12, 1, false, false, false, 'No response big SNR';...
    'desmond16', '20200312', 8, 1, false, false, false, 'No response big SNR';...
    'desmond16', '20200312', 3, 1, false, false, false, 'No response big SNR';...
    'desmond16', '20200312', 1, 1, false, false, false, 'No response big SNR';...
    'desmond16', '20200314', 1, 1, false, false, false, 'No response big SNR';...
    'desmond16', '20200314', 2, 1, false, false, false, 'No response big SNR';...
    'desmond16', '20200314', 3, 1, false, false, false, 'No response big SNR';...
    'desmond16', '20200314', 6, 1, false, false, true, 'Complex response big SNR no collision';...
    'desmond16', '20200314', 7, 1, false, false, true, 'Complex response big SNR no collision';...
    'desmond16', '20200314', 10, 1, false, false, true, 'Complex response big SNR no collision';...
    'desmond16', '20200314', 11, 1, false, false, true, 'Complex response big SNR no collision';...
    'desmond16', '20200314', 12, 1, false, false, true, 'Complex response big SNR no collision';...
    'desmond16', '20200315', 6, 1, false, false, true, 'Complex response big SNR no collision';...
    'desmond16', '20200315', 4, 1, false, false, true, 'Complex response big SNR no collision';...
    'desmond16', '20200316', 16, 1, false, false, true, 'Complex response big SNR no collision';...
    'desmond16', '20200316', 12, 1, false, false, false, 'No response big SNR no collision';...
    'desmond16', '20200316', 10, 1, false, false, false, 'No response big SNR no collision';...
    'desmond16', '20200316', 9, 1, false, false, true, 'Complex response big SNR no collision';...
    'desmond16', '20200316', 8, 1, false, false, false, 'No response med SNR no collision';...
    'desmond16', '20200316', 7, 1, false, false, false, 'No response small SNR no collision';...
    'desmond16', '20200316', 1, 1, false, false, false, 'No response small SNR no collision';...
    'desmond16', '20200317', 1, 1, false, false, false, 'No response small SNR no collision';...
    'desmond16', '20200317', 2, 1, false, false, false, 'No response small SNR no collision';...
    'desmond16', '20200317', 5, 1, false, false, false, 'No response small SNR no collision';...
    'desmond17', '20200225', 2, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200225', 4, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200225', 6, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200225', 7, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200225', 10, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200225', 11, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200225', 13, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200227', 5, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200227', 4, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200227', 3, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200227', 1, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200228', 1, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200228', 3, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200228', 4, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200228', 5, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200228', 6, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200228', 7, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200228', 8, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200302', 1, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200302', 2, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200302', 6, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200302', 7, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200302', 8, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200303', 1, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200303', 2, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200303', 4, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200303', 9, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200303', 10, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200304', 1, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200304', 3, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200304', 9, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200304', 11, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200304', 13, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200304', 14, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200305', 2, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200305', 3, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200305', 4, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200305', 5, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200305', 6, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200305', 7, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200306', 1, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200306', 2, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200306', 3, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200306', 4, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200306', 5, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200306', 7, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200306', 10, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200306', 11, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200306', 12, 1, false, true, true, '2ms response, but no collision';...
    'desmond17', '20200309', 1, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200309', 2, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200309', 4, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200309', 5, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200309', 3, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200309', 7, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200309', 8, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200309', 9, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200310', 2, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200310', 3, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200310', 5, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200310', 8, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200310', 10, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200311', 1, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200311', 2, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200311', 3, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200311', 4, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200311', 5, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200311', 6, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200312', 1, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200312', 6, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200315', 1, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200315', 4, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200315', 5, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200315', 7, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200315', 8, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200316', 2, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200316', 5, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200316', 8, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200316', 9, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200316', 10, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200316', 11, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200317', 1, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200317', 5, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200317', 8, 1, false, false, false, 'No response no collision';...
    'desmond17', '20200317', 10, 1, false, false, false, 'No response no collision. Desmond17 seems to have plenty of big units, but none had stim response.';...
    'desmond18', '20200224', 7, 1, false, false, false, 'Maybe a 10ms stim response? .';...
    'desmond18', '20200224', 2, 1, false, false, false, 'No response.';...
    'desmond18', '20200224', 4, 1, false, false, false, 'No response.';...
    'desmond18', '20200224', 5, 1, false, false, false, 'No response.';...
    'desmond18', '20200224', 6, 1, false, false, false, 'No response.';...
    'desmond18', '20200224', 7, 1, false, false, false, 'No response.';...
    'desmond18', '20200224', 8, 1, false, false, false, 'No response.';...
    'desmond18', '20200224', 9, 1, false, false, false, 'No response.';...
    'desmond18', '20200224', 10, 1, false, false, false, 'No response.';...
    'desmond18', '20200224', 11, 1, false, false, false, 'No response.';...
    'desmond18', '20200224', 12, 1, false, false, false, 'No response.';...
    'desmond18', '20200224', 13, 1, false, false, false, 'No response.';...
    'desmond18', '20200224', 15, 1, false, false, false, 'No response.';...
    'desmond18', '20200225', 16, 1, false, false, true, 'Complex Stim response no collision high SNR.';...
    'desmond18', '20200225', 10, 1, false, false, true, 'Complex Stim response no collision high SNR.';...
    'desmond18', '20200225', 9, 1, false, false, true, 'Complex Stim response no collision high SNR.';...
    'desmond18', '20200225', 8, 1, false, false, true, 'Complex Stim response no collision high SNR.';...
    'desmond18', '20200225', 4, 1, false, false, true, 'Complex Stim response no collision high SNR.';...
    'desmond18', '20200225', 2, 1, false, false, true, 'Complex Stim response no collision high SNR.';...
    'desmond18', '20200226', 1, 1, false, false, true, 'Complex Stim response no collision high SNR.';...
    'desmond18', '20200226', 3, 1, false, false, true, 'Complex Stim response no collision high SNR.';...
    'desmond18', '20200226', 4, 1, false, false, true, 'Complex Stim response no collision high SNR.';...
    'desmond18', '20200227', 1, 1, true, true, true, 'Stim response, collisions. Low SNR.';...
    'desmond18', '20200227', 2, 1, false, false, true, 'Complex Stim response no collision high SNR.';...
    'desmond18', '20200227', 3, 1, false, false, true, 'Complex Stim response no collision high SNR.';...
    'desmond18', '20200227', 7, 1, false, false, true, 'Complex Stim response no collision high SNR.';...
    'desmond18', '20200227', 8, 1, false, false, true, 'Complex Stim response no collision high SNR.';...
    'desmond18', '20200227', 9, 1, false, false, true, 'Complex Stim response no collision high SNR.';...
    'desmond18', '20200227', 12, 1, false, false, true, 'Complex Stim response no collision high SNR.';...
    'desmond18', '20200227', 13, 1, false, false, false, 'No Stim response no collision high SNR.';...
    'desmond18', '20200227', 14, 1, false, false, false, 'no stim response no collision high SNR.';...
    'desmond18', '20200228', 1, 1, false, false, false, 'No stim response no collision high SNR.';...
    'desmond18', '20200228', 2, 1, false, false, false, 'No stim response no collision high SNR.';...
    'desmond18', '20200228', 7, 1, false, false, false, 'No stim response no collision high SNR.';...
    'desmond18', '20200228', 8, 1, false, false, false, 'No Stim response no collision.';...
    'desmond18', '20200228', 11, 1, false, false, false, 'no stim response no collision.';...
    'desmond18', '20200302', 1, 1, false, false, false, 'No stim response no collision high SNR.';...
    'desmond18', '20200302', 6, 1, false, false, false, 'No stim response no collision high SNR.';...
    'desmond18', '20200302', 4, 1, false, false, false, 'No stim response no collision low SNR.';...
    'desmond18', '20200302', 7, 1, false, false, false, 'No stim response no collision high SNR.';...
    'desmond18', '20200303', 2, 1, false, false, false, 'No stim response no collision low SNR.';...
    'desmond18', '20200303', 5, 1, false, false, false, 'No stim response no collision high SNR.';...
    'desmond18', '20200303', 8, 1, false, false, false, 'No stim response no collision low SNR.';...
    'desmond18', '20200303', 1, 1, false, false, false, 'No stim response no collision high SNR.';...
    'desmond18', '20200303', 3, 1, false, false, false, 'No stim response no collision low SNR. TODO: Check 2ms stims later.';...
    'desmond18', '20200303', 5, 1, false, false, false, 'No stim response no collision high SNR.';...
    'desmond18', '20200303', 6, 1, false, false, false, 'No stim response no collision low SNR.';...
    'desmond18', '20200303', 7, 1, false, false, false, 'No stim response no collision high SNR.';...
    'desmond18', '20200303', 10, 1, false, false, false, 'No stim response no collision low SNR.';...
    'desmond18', '20200303', 12, 1, false, false, false, 'No stim response no collision high SNR.';...
    'desmond18', '20200306', 1, 1, true, false, true, 'Stim response, maybe collision but all long trains high SNR.';...
    'desmond18', '20200306', 2, 1, false, false, false, 'No stim response no collision high SNR.';...
    'desmond18', '20200306', 4, 1, false, false, false, 'No stim response no collision high SNR.';...
    'desmond18', '20200306', 5, 1, false, false, true, 'Stim response, maybe collision but all long trains high SNR.';...
    'desmond18', '20200306', 6, 1, true, false, true, 'Stim response, maybe collision but all long trains high SNR.';...
    'desmond18', '20200306', 7, 1, false, false, false, 'No stim response no collision high SNR.';...
    'desmond18', '20200306', 9, 1, true, false, true, 'Stim response, maybe collision but all long trains high SNR.';...
    'desmond18', '20200310', 1, 1, false, true, true, 'Stim response, no collision high SNR.';...
    'desmond18', '20200310', 2, 1, false, false, true, 'Complex Stim response, no collisionn high SNR.';...
    'desmond18', '20200310', 6, 1, false, false, true, 'Complex Stim response, no collision low SNR.';...
    'desmond18', '20200310', 7, 1, false, false, true, 'Stim response but big temporal jitter, high SNR.';...
    'desmond18', '20200310', 8, 1, false, false, true, 'Stim response, maybe collision low SNR.';...
    'desmond18', '20200311', 1, 1, false, false, false, 'No Stim response, no collision high SNR.';...
    'desmond18', '20200311', 2, 1, false, false, true, 'Stim response but big temporal jitter, high SNR.';...
    'desmond18', '20200311', 5, 1, true, true, true, 'Stim response, Collisions. high SNR..';...
    'desmond18', '20200311', 6, 1, false, false, true, 'Complex Stim response, no collision high SNR.';...
    'desmond18', '20200311', 8, 1, false, false, false, 'No Stim response, no collision high SNR.';...
    'desmond18', '20200311', 9, 1, false, true, true, 'Stim response, no collision high SNR.';...
    'desmond18', '20200311', 10, 1, false, false, true, 'Stim response but big temporal jitter, high SNR.';...
    'desmond18', '20200314', 1, 1, false, true, true, 'Stim response, no collisions, high SNR';...
    'desmond18', '20200314', 6, 1, false, false, false, 'No stim response, no collisions, high SNR';...
    'desmond18', '20200314', 7, 1, false, false, false, 'No stim response, no collisions, high SNR';...
    'desmond18', '20200314', 9, 1, false, false, false, 'No stim response, no collisions, high SNR';...
    'desmond18', '20200314', 10, 1, false, true, true, 'Stim response, no collisions, high SNR';...
    'desmond18', '20200315', 1, 1, false, false, true, 'Complex stim response, no collisions, high SNR';...
    'desmond18', '20200315', 2, 1, false, false, true, 'Complex stim response, no collisions, low SNR';...
    'desmond18', '20200315', 3, 1, false, false, true, 'Complex stim response, no collisions, low SNR';...
    'desmond18', '20200315', 6, 1, false, false, false, 'No Stim response, no collisions, high SNR';...
    'desmond18', '20200315', 7, 1, false, false, false, 'No Stim response, no collisions, high SNR';...
    'desmond18', '20200315', 8, 1, false, false, false, 'No Stim response, no collisions, high SNR';...
    'desmond18', '20200315', 9, 1, true, true, true, 'Stim response, Collisions, high SNR';...
    'desmond18', '20200315', 10, 1, false, false, false, 'No stim response, no collisions, high SNR';...
    'desmond18', '20200315', 11, 1, false, false, false, 'No Stim response, no collisions, high SNR'...
    };

    ctUnitDesc = struct('Animal', collisionTestUnits(:, 1), 'Date', collisionTestUnits(:, 2), 'Channel', collisionTestUnits(:, 3), 'Unit', collisionTestUnits(:, 4), 'HasCollision', collisionTestUnits(:, 5), 'HasTrueResponse', collisionTestUnits(:, 6), 'HasAnyResponse', collisionTestUnits(:, 7), 'Notes', collisionTestUnits(:, 8));

ctud = CollisionTestUnitDesc(collisionTestUnits);
ctud = ctud([ctud.HasAnyResponse]);
% ctud.read();

selCollisionTestResults = {...


};





% Generate PETH data struct
expNames = cell(length(batchPlotList), 1);
for iExp = 1:length(batchPlotList)
    expNames{iExp} = [batchPlotList{iExp, 1}, '_', num2str(batchPlotList{iExp, 2})];
end

expNamesUnique = unique(expNames);

for iTr = 1:length(expNamesUnique)
    tr = TetrodeRecording.BatchLoad(expNamesUnique(iTr));
    try
        if iTr == 1;
            PETH = TetrodeRecording.BatchPETHistCounts(tr, batchPlotList, 'TrialLength', 6, 'ExtendedWindow', 1, 'SpikeRateWindow', 100, 'ExtendedWindowStim', [-1, 1], 'SpikeRateWindowStim', 10, 'Press', true, 'Lick', false, 'Stim', true);
        else
            PETH = [PETH, TetrodeRecording.BatchPETHistCounts(tr, batchPlotList, 'TrialLength', 6, 'ExtendedWindow', 1, 'SpikeRateWindow', 100, 'ExtendedWindowStim', [-1, 1], 'SpikeRateWindowStim', 10, 'Press', true, 'Lick', false, 'Stim', true)];
        end
    catch ME
        warning(['Error when processing iTr = ', num2str(iTr), ' - this one will be skipped.'])
        warning(sprintf('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message))
    end
end