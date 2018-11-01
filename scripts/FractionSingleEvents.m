function FractionSingleEventsOutput = FractionSingleEvents(currCell)

    events = EventTimes(currCell);
    nEvents = numel(events);
    
    burstEvents = EventsPerBurst(currCell);
    nBurstEvents = sum(burstEvents.nEvents);
    
    FractionSingleEventsOutput = (nEvents - nBurstEvents)/nEvents;

end