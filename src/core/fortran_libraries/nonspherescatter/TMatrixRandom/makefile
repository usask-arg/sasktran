#-----------------------------------------------------------------------------
# Thse macros are normally defined by NMAKE on the command line
# OBJDIR  = .\obj	Location to write obj files
# SRCDIR  = .\      Location to find source files
# LIBDIR  = .\      Location to write *.lib and *.dll files
#-----------------------------------------------------------------------------

OBJDIR  = ..\..\obj\tmatrix
#SRCDIR  = .\
LIBDIR  = ..\..\repos_libs\x64

FC      = lf95
FFLAGS  = -winconsole -c -wide -w -nco -dbl -o1 -ml msvc -MOD:$(OBJDIR) -tpp -chk(a,s)
LFLAGS  = -nco -nomap -nologo 

.SUFFIXES:
.SUFFIXES: .obj .f
{$(SRCDIR)}.f{$(OBJDIR)}.obj:
	$(FC) $(FFLAGS) -o $@ $**
	
$(LIBDIR)\tmatrix.dll: $(OBJDIR)\tmd.lp.obj $(OBJDIR)\lpd.obj
    attrib -r $(LIBDIR)\tmatrix.*
	$(FC) $(LFLAGS) -dll -out $@ $**

#test: $(OBJDIR)\MVTstNew.obj $(LIBDIR)\tmatrix.dll $(OBJDIR)\lpd.obj
#	$(FC) $(LFLAGS) -out MVTstNew.exe $** $(LIBDIR)\tmatrix.lib
	
clean:
#    attrib -r $(LIBDIR)\tmatrix.*
#	del /q /f $(OBJDIR)\MVTstNew.obj
#	del /q /f $(OBJDIR)\tmatrix.obj
