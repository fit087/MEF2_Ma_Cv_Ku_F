# Microsoft Developer Studio Project File - Name="trie2d" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=trie2d - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "trie2d.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "trie2d.mak" CFG="trie2d - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "trie2d - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "trie2d - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "trie2d - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD BASE RSC /l 0x416 /d "NDEBUG"
# ADD RSC /l 0x416 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "trie2d - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# SUBTRACT F90 /check:bounds
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD BASE RSC /l 0x416 /d "_DEBUG"
# ADD RSC /l 0x416 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib /nologo /subsystem:console /incremental:no /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "trie2d - Win32 Release"
# Name "trie2d - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Source File

SOURCE=..\ADDBAN2.F
DEP_F90_ADDBA=\
	"..\tapes.h"\
	
# End Source File
# Begin Source File

SOURCE=..\ASSEMB.F
DEP_F90_ASSEM=\
	"..\common.h"\
	
# End Source File
# Begin Source File

SOURCE=..\BKDATA.F
DEP_F90_BKDAT=\
	"..\cntl.h"\
	
# End Source File
# Begin Source File

SOURCE=..\COLSOL.F
# End Source File
# Begin Source File

SOURCE=..\CONTRL.F
DEP_F90_CONTR=\
	"..\cntl.h"\
	"..\tapes.h"\
	
# End Source File
# Begin Source File

SOURCE=..\EQUAT.F
DEP_F90_EQUAT=\
	"..\cntl.h"\
	"..\tapes.h"\
	
# End Source File
# Begin Source File

SOURCE=..\FUNCOES.F
# End Source File
# Begin Source File

SOURCE=..\INELEM.F
DEP_F90_INELE=\
	"..\tapes.h"\
	
# End Source File
# Begin Source File

SOURCE=..\INMESH.F
DEP_F90_INMES=\
	"..\cntl.h"\
	"..\common.h"\
	"..\pointers.h"\
	
# End Source File
# Begin Source File

SOURCE=..\IOMNGR.F
DEP_F90_IOMNG=\
	"..\tapes.h"\
	
# End Source File
# Begin Source File

SOURCE=..\JOINTS.F
DEP_F90_JOINT=\
	"..\cntl.h"\
	"..\tapes.h"\
	
# End Source File
# Begin Source File

SOURCE=..\LDRIVER.F
DEP_F90_LDRIV=\
	"..\cntl.h"\
	"..\common.h"\
	"..\pointers.h"\
	
# End Source File
# Begin Source File

SOURCE=..\LOADS.F
DEP_F90_LOADS=\
	"..\cntl.h"\
	"..\common.h"\
	"..\tapes.h"\
	
# End Source File
# Begin Source File

SOURCE=..\MATRIZ.F
DEP_F90_MATRI=\
	"..\cntl.h"\
	"..\common.h"\
	"..\pointers.h"\
	"..\tapes.h"\
	
# End Source File
# Begin Source File

SOURCE=..\MATSET.F
DEP_F90_MATSE=\
	"..\tapes.h"\
	
# End Source File
# Begin Source File

SOURCE=..\MEMO.F
DEP_F90_MEMO_=\
	"..\common.h"\
	
# End Source File
# Begin Source File

SOURCE=..\PROFIL2.F
DEP_F90_PROFI=\
	"..\tapes.h"\
	
# End Source File
# Begin Source File

SOURCE=..\PTLOAD.F
DEP_F90_PTLOA=\
	"..\tapes.h"\
	
# End Source File
# Begin Source File

SOURCE=..\QUAD4.f
DEP_F90_QUAD4=\
	"..\tapes.h"\
	
# End Source File
# Begin Source File

SOURCE=..\SAIDAS.F
DEP_F90_SAIDA=\
	"..\cntl.h"\
	"..\tapes.h"\
	
# End Source File
# Begin Source File

SOURCE=..\SETDOF.F
# End Source File
# Begin Source File

SOURCE=..\SOLVER.F
DEP_F90_SOLVE=\
	"..\cntl.h"\
	"..\common.h"\
	"..\pointers.h"\
	"..\tapes.h"\
	
# End Source File
# Begin Source File

SOURCE=..\STRESS.F
DEP_F90_STRES=\
	"..\tapes.h"\
	
# End Source File
# Begin Source File

SOURCE=..\TIMETB.F
DEP_F90_TIMET=\
	"..\tapes.h"\
	
# End Source File
# Begin Source File

SOURCE=..\TRIE2D2.F
DEP_F90_TRIE2=\
	"..\tapes.h"\
	
# End Source File
# Begin Source File

SOURCE=..\TRINL.F
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
