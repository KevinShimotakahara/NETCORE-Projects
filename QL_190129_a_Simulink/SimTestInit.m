global UDs
sn_FieldLength = 5;
sizePDU = 772;
UDs(1).RLCtxEntity = RLCentityUMtx(sn_FieldLength,sizePDU);
UDs(5).RLCtxEntity = RLCentityUMtx(sn_FieldLength,sizePDU);
UDs(7).RLCtxEntity = RLCentityUMtx(sn_FieldLength,sizePDU);

UDs(2).RLCrxEntity = RLCentityUMrx(sn_FieldLength,sizePDU);
UDs(3).RLCrxEntity = RLCentityUMrx(sn_FieldLength,sizePDU);
UDs(4).RLCrxEntity = RLCentityUMrx(sn_FieldLength,sizePDU);
UDs(5).RLCrxEntity = RLCentityUMrx(sn_FieldLength,sizePDU);
UDs(6).RLCrxEntity = RLCentityUMrx(sn_FieldLength,sizePDU);

for k = 1:7
    UDs(k).nodeID = k;
    UDs(k).pktNumber = [];
end