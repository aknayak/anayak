#! /bin/csh
foreach mass (80 100 120 140 150 155 160)

echo  "processing mH${mass}"

../test/lands.exe --PhysicsModel ChargedHiggs -M Hybrid --bQuickEstimateInitialLimit 0 --initialRmin 0. --initialRmax 0.2 -d datacard_csbar_mH${mass}.txt  --doExpectation 1 -t 100 > & datacard_csbar_mH${mass}.txt_out.txt

end
