function TimeFirstBurstOutput = TimeFirstBurst(currCell)

firstBurstTime = currCell.Bursts(1,2);
TimeFirstBurstOutput = firstBurstTime - currCell.startTime;



end