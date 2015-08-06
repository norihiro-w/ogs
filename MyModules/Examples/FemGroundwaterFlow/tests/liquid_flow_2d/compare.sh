echo "# ogs6-thmc"
time /home/localadmin/ogs/ogs6/ogs6-thmc/BuildClangRelease/bin/ogs6 q_quad -i 2 -p 0 -initx_zeros 1 &> ogs.log
#echo ""
#echo "# ogs6-eigen"
#time /home/localadmin/ogs/ogs6/ogs6-thmc/BuildClangRelease/bin/ogs6-eigen q_quad -i 2 -p 0 -initx_zeros 1 &> ogs.log
echo ""
echo "#ogs5-mkl"
time /home/localadmin/ogs/ogs5/trunk/trunk-sources-git-ufz/BuildMKLRelease/bin/ogs q_quad &> ogs5.log
#time ~/Work/ogs/ogs6/ogs6-thmc/BuildClangRelease/bin/ogs6 q_quad &> ogs.log
#time ~/Work/ogs/ogs5/trunk-github/BuildGnuMKLRelease/bin/ogs q_quad &> ogs5.log
