# installing/loading the package:
if(!require(installr)) { install.packages("installr"); require(installr)} #load / install+load installr

updateR(F, T, T, F, T, F, T) # install, move, update.package, quit R.
updateR(T, T, T, T, T, T, T)
# # the safest upgrade option: See the NEWS,
# install R, copy packages, keep old packages,
# update packages in the new installation,
# start the Rgui of the new R, and quite current session
# of R
updateR()#will ask at every stage of the process



# to see if there are any issues with latest version
inst <- packageStatus()$inst
inst[inst$Status != "ok", c("Package", "Version", "Status")]