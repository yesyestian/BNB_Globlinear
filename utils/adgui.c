/*
   File:            spr.c
   Author:          Daniel Halstead
   Created:         26th Apr 1999
   Modified:        

   Description:     Andrews Graphics
*/

#include "standard.h"
#include "adgui.h"
#ifdef PC_MVIS_PLATFORM
#include "amgr.h"
#include "amxw.h"
#endif
#include "amiv.h"

FILE *LogFile;                  /*Logs all output to the text window, if not NULL*/
int IShallDie;                  /*Indicates to terminate current loop*/

#ifdef PC_MVIS_PLATFORM

#include <process.h>

/*Here we need to declare a lot of global variables.  This is because
  these are variables that both threads will need to access.*/

int Processing; /*Holds the state of the input window*/
HANDLE PressedEnter,SwitchedStatus,KeyPress,MouseClick,ClickOrEnter,FinishedCommand; /*Events*/
CRITICAL_SECTION critsect;
HWND hWnd,hEditIn,hEditOut,hPict,hEditProc;                /*Windows*/
HINSTANCE hInstance;
OPENFILENAME OFN;               /*Holds parameters for file dialogs*/
int NumClientLines;             /*The number of text-lines in the window*/
int clientwidth,clientheight;   /*Dimensions of window in pixels*/
int panewidth,paneheight;       /*Dimensions of graphic pane in logical units*/
int MemWidth=512,MemHeight=512; /*Dimensions of graphic being stored in memory*/
int xScroll=0,yScroll=0;        /*Top-left coordinate of graphics as it is being displayed*/
bool ScrollbarsOn;              /*Whether scrollbars are in graphic panel*/
int charwidth,charheight;       /*Dimensions of a character in pixels*/
int GFontSize = 2;              /*The size of text in the left pane.*/
int TFontSize = 1;              /*The size of text in the right pane.*/
int XSplit;                     /*Where the split occurs between EditIn and EditProc*/
char OutBuf[BUFSIZE];           /*The buffer for the text window*/
char EditInBuf[LINESIZE];       /*The buffer for the input window*/
char GetcBuf;                   /*Holds the last key pressed.*/
int FlushFlag;                  /*Indicates to flush output*/
POINT MouseCoords;              /*Holds the coords of the last click*/
int MouseButton;                /*Holds the last button pressed: 1=left, 3=right*/
void *DlgInfo;                  /*Contains any info the programmer wants a dialog to know about*/
int main_argc;                  /*The number of arguments to main()*/
char **main_argv;               /*The arguments going to main()*/

#define ID_EditIn 1
#define ID_EditOut 2
#define ID_Pict 3
#define ID_EditProc 4
#define DEF_PICT_STYLE WS_CHILD|WS_VISIBLE|SS_WHITEFRAME

extern void main(int argc,char *argv[]);

extern HDC hPictDC;
extern HBITMAP hBitmap;
extern void old_screen_off(void);

void pre_main(void* pvoid);
void init_OFN(HWND hwnd);
int screenprint(char *outstr);
void UpdateEditOut();
void gui_run_menu(int id);

LRESULT CALLBACK WndProc (HWND, UINT, WPARAM, LPARAM) ;
LRESULT CALLBACK PictWndProc (HWND, UINT, WPARAM, LPARAM) ;
LRESULT CALLBACK EditProcProc (HWND, UINT, WPARAM, LPARAM) ;
LRESULT CALLBACK EditInProc (HWND, UINT, WPARAM, LPARAM) ;
LRESULT CALLBACK PictProc (HWND, UINT, WPARAM, LPARAM) ;
WNDPROC OldEditInProc;
WNDPROC OldEditProcProc;
WNDPROC OldPictProc;

typedef struct menu_item{
  char name[50];
  int id;
} menu_item;

typedef struct pict_info{
  int height,width;
  HDC hmemdc;
  HWND hwnd;
  HBITMAP hbitmap;
  struct pict_info *next;
} pict_info;

#ifdef DO_MENU
menu_item MenuCommand[MAX_MENU-MIN_MENU+1];
#else
menu_item *MenuCommand = NULL;
#endif
void (*Dlg_Dispatch)(void *env,void *c,bool act) = NULL;
pict_info *pict_list = NULL;
int curr_command = 0,last_command = -1;
char command_hist[COMHISTSIZE][LINESIZE];
int CaretPos = 0, NextLine = 0, NumOutLines = 0;

const char szAppName[] = "SPR Console" , szPictName[] = "SPR Graphic";

WPARAM message_loop(void* pvoid){
  MSG msg ;
  while (GetMessage (&msg, NULL, 0, 0)){
    TranslateMessage (&msg) ;
	  DispatchMessage (&msg) ;
  }
  return msg.wParam;
}

int WINAPI WinMain (HINSTANCE hInst, HINSTANCE hPrevInstance,
                    PSTR szCmdLine, int iCmdShow){
    WNDCLASSEX  wndclass, pictwndclass;
    string_array *sa = NULL;
    int i,maxlen = 0;
    int window_height;
    char caption[50];
    hInstance = hInst;

    /*Define two window classes: the main window, and the pop-up picture windows*/
    wndclass.cbSize        = sizeof (wndclass) ;
    wndclass.style         = CS_HREDRAW | CS_VREDRAW;
    wndclass.lpfnWndProc   = WndProc ;
    wndclass.cbClsExtra    = 0 ;
    wndclass.cbWndExtra    = 0 ;
    wndclass.hInstance     = hInstance ;
    wndclass.hIcon         = LoadIcon (hInstance, MAKEINTRESOURCE(IDI_SPR)) ;
    wndclass.hCursor       = LoadCursor (NULL, IDC_ARROW) ;
    wndclass.hbrBackground = (HBRUSH) GetStockObject (WHITE_BRUSH) ;
    wndclass.lpszClassName = szAppName ;
    wndclass.hIconSm       = LoadIcon (NULL, NULL) ;
#ifdef DO_MENU
    wndclass.lpszMenuName  = MAKEINTRESOURCE(IDR_MENU) ;
#else
    wndclass.lpszMenuName  = NULL ;
#endif

    pictwndclass.cbSize        = sizeof (pictwndclass) ;
    pictwndclass.style         = CS_HREDRAW | CS_VREDRAW;
    pictwndclass.lpfnWndProc   = PictWndProc ;
    pictwndclass.cbClsExtra    = 0 ;
    pictwndclass.cbWndExtra    = 0 ;
    pictwndclass.hInstance     = hInstance ;
    pictwndclass.hIcon         = LoadIcon (NULL, IDI_APPLICATION) ;
    pictwndclass.hCursor       = LoadCursor (NULL, IDC_ARROW) ;
    pictwndclass.hbrBackground = (HBRUSH) GetStockObject (WHITE_BRUSH) ;
    pictwndclass.lpszMenuName  = NULL ;
    pictwndclass.lpszClassName = szPictName ;
    pictwndclass.hIconSm       = LoadIcon (NULL, IDI_APPLICATION) ;

    /*Set up the arguments to be passed into main()*/
    sa = mk_broken_string(szCmdLine);
    insert_in_string_array(sa,0,"adgui.exe");
    main_argc = string_array_size(sa);
    if(main_argc>=3 && caseless_eq_string(string_array_ref(sa,1),"gfont")){
      GFontSize = atoi(string_array_ref(sa,2));
      main_argc-=2;
      string_array_remove(sa,2);
      string_array_remove(sa,1);
    }
    if(main_argc>=3 && caseless_eq_string(string_array_ref(sa,1),"tfont")){
      TFontSize = atoi(string_array_ref(sa,2));
      main_argc-=2;
      string_array_remove(sa,2);
      string_array_remove(sa,1);
    }
    main_argv = (char**) malloc(sizeof(char*)*main_argc);
    for(i=0;i<main_argc;i++){
      main_argv[i] = (char*) malloc(sizeof(char)*(strlen(string_array_ref(sa,i))+1));
      strcpy(main_argv[i],string_array_ref(sa,i));
    } 

    /*Tell the registry about our new window classes*/
    RegisterClassEx (&wndclass) ;
    RegisterClassEx (&pictwndclass) ;

    /*Initialize some globals*/
    IShallDie = FALSE;
    XSplit = -1;
    LogFile = NULL;
#ifdef DO_MENU
    for(i=0;i<=(MAX_MENU-MIN_MENU);i++)
      MenuCommand[i].id = 0;
    window_height = 464;           // initial y size    -- could be 556
#else
    window_height = 444;           // initial y size    -- could be 556
#endif

    /*Create the main window based on this class, and initialize it to be busy processing*/
    Processing = CLI_PROCESSING;
    LoadString(NULL,IDS_CAPTION,caption,50);
    hWnd = CreateWindow (szAppName,         // window class name
						caption,     // window caption
						WS_OVERLAPPEDWINDOW,     // window style
	          CW_USEDEFAULT,           // initial x position
		        CW_USEDEFAULT,           // initial y position
			      800,           // initial x size  -- could be 1024
            window_height,
					  NULL,                    // parent window handle
						NULL,                    // window menu handle
	          hInstance,               // program instance handle
			      NULL) ;		             // creation parameters
    ShowWindow (hWnd, iCmdShow) ;
    UpdateWindow (hWnd) ;
    free_string_array(sa);
    init_OFN(hWnd);

    /*Initialize the command history and the stdout flushing*/
    for(i=0;i<COMHISTSIZE;i++) strcpy(command_hist[i],"");
    FlushFlag = FALSE;

    /*The big daddy message loop*/
    return message_loop(NULL);

    //_cexit();
    //abort();
}

//#endif

void pre_main(void* pvoid){
  main(main_argc,main_argv);
  gui_quit();
}

LRESULT CALLBACK WndProc (HWND hwnd, UINT iMsg, WPARAM wParam, LPARAM lParam){
	static HBRUSH hBrushDeactive;
	static HFONT hFontDeactive;
  LOGFONT lf;
  HDC hdc ;
  PAINTSTRUCT ps ;
	TEXTMETRIC tm;
  SCROLLINFO si;
	int id;

  switch (iMsg){
	case WM_CREATE :
		hdc = GetDC(hwnd);
		GetTextMetrics(hdc,&tm);
		charwidth = tm.tmAveCharWidth;
		charheight = tm.tmHeight+tm.tmExternalLeading;
		ReleaseDC(hwnd,hdc);

//		InitCommonControls();
		hEditProc = CreateWindowEx(WS_EX_OVERLAPPEDWINDOW,"edit",NULL,
							WS_CHILD|WS_VISIBLE|ES_LEFT|ES_READONLY,
							0,0,0,0,hwnd,(HMENU) ID_EditProc,
							((LPCREATESTRUCT) lParam)->hInstance,NULL);
		hEditIn = CreateWindowEx(WS_EX_OVERLAPPEDWINDOW,"edit",NULL,
							WS_CHILD|ES_LEFT|ES_AUTOHSCROLL|ES_WANTRETURN,
							0,0,0,0,hwnd,(HMENU) ID_EditIn,
							((LPCREATESTRUCT) lParam)->hInstance,NULL);
		hEditOut = CreateWindowEx(WS_EX_OVERLAPPEDWINDOW,"edit",NULL,
							WS_CHILD|WS_VISIBLE|WS_HSCROLL|WS_VSCROLL|
							ES_LEFT|ES_MULTILINE|ES_AUTOVSCROLL|ES_READONLY,
							0,0,0,0,hwnd,(HMENU) ID_EditOut,
							((LPCREATESTRUCT) lParam)->hInstance,NULL);
		hPict = CreateWindowEx(WS_EX_STATICEDGE,"static",NULL,
							WS_CHILD|WS_VISIBLE|SS_WHITEFRAME|WS_HSCROLL|WS_VSCROLL,
							0,0,0,0,hwnd,(HMENU) ID_Pict,
							((LPCREATESTRUCT) lParam)->hInstance,NULL);

		OldEditProcProc = (WNDPROC) SetWindowLong(hEditProc,GWL_WNDPROC,(LONG) EditProcProc);
		OldEditInProc = (WNDPROC) SetWindowLong(hEditIn,GWL_WNDPROC,(LONG) EditInProc);
		OldPictProc = (WNDPROC) SetWindowLong(hPict,GWL_WNDPROC,(LONG) PictProc);

		hBrushDeactive = CreateSolidBrush(GetSysColor(COLOR_BTNFACE));
		hFontDeactive = (HFONT) GetStockObject(OEM_FIXED_FONT);
    GetObject(hFontDeactive,sizeof(LOGFONT),&lf);
    lf.lfHeight = lf.lfHeight-2+TFontSize;
    lf.lfWidth = lf.lfWidth-3+TFontSize;
    hFontDeactive = CreateFontIndirect(&lf);
		SendMessage(hEditOut,WM_SETFONT,(WPARAM) hFontDeactive,FALSE);

    si.nPage = 512;
    si.nMin = 0;
    si.nMax = 1000;
    si.fMask = SIF_PAGE|SIF_RANGE;
    SetScrollInfo(hPict,SB_HORZ,&si,FALSE);
    SetScrollInfo(hPict,SB_VERT,&si,FALSE);
    ShowScrollBar(hPict,SB_BOTH,FALSE);

    SetWindowText(hEditProc,"Processing....");
    SetFocus(hEditIn);

    PressedEnter = CreateEvent(NULL,FALSE,FALSE,NULL);
    SwitchedStatus = CreateEvent(NULL,FALSE,FALSE,NULL);
    KeyPress = CreateEvent(NULL,FALSE,FALSE,NULL);
    MouseClick = CreateEvent(NULL,FALSE,FALSE,NULL);
    ClickOrEnter = CreateEvent(NULL,FALSE,FALSE,NULL);
    FinishedCommand = CreateEvent(NULL,FALSE,FALSE,NULL);
    InitializeCriticalSection(&critsect);

    _beginthread(pre_main,0,NULL);
		return 0 ;

	case WM_SIZE:
		clientwidth = LOWORD(lParam);
		clientheight = HIWORD(lParam);
    if(XSplit==-1) XSplit = clientwidth;
		NumClientLines = (int) (clientheight-charheight*1.5)/(charheight-2);

		hdc = GetDC(hwnd);
    MoveWindow(hEditProc,0,(int) (clientheight-charheight*1.5),XSplit,(int) (charheight*1.5),TRUE);
    MoveWindow(hEditIn,XSplit,(int) (clientheight-charheight*1.5),clientwidth,(int) (charheight*1.5),TRUE);
		MoveWindow(hEditOut,clientwidth/2,0,clientwidth/2,(int) (clientheight-charheight*1.5),TRUE);
		MoveWindow(hPict,0,0,clientwidth/2,(int) (clientheight-charheight*1.5),TRUE);
		/*MoveWindow(hPict,0,0,clientwidth/2-GetSystemMetrics(SM_CXVSCROLL),
                    (int) (clientheight-charheight*1.5-GetSystemMetrics(SM_CYHSCROLL)),TRUE);*/
    if(ScrollbarsOn){
      panewidth = clientwidth/2-2-GetSystemMetrics(SM_CXVSCROLL);
      paneheight = (int) (clientheight-charheight*1.5-GetSystemMetrics(SM_CYHSCROLL));
    } else {
      panewidth = clientwidth/2-2;
      paneheight = (int) (clientheight-charheight*1.5);
    }
    if(panewidth==paneheight) panewidth = paneheight = 512;
    else if(paneheight>panewidth){
      paneheight = 512 * paneheight / panewidth;
      panewidth = 512;
    } else {
      panewidth = 512 * panewidth / paneheight;
      paneheight = 512;
    }
    if(ScrollbarsOn) SendMessage(hPict,WM_SHOWSCROLL,0,ScrollbarsOn);

		ReleaseDC(hwnd,hdc);
		break;

	case WM_PAINT :
		hdc = BeginPaint (hwnd, &ps) ;
		RedrawPict(hdc);
		EndPaint (hwnd, &ps) ;
    return 0 ;

	case WM_CTLCOLOREDIT :
		id = GetWindowLong((HWND) lParam,GWL_ID);
		if(id==ID_EditOut){
			SetBkColor((HDC) wParam,GetSysColor(COLOR_BTNFACE));
			SetTextColor((HDC) wParam,GetSysColor(COLOR_BTNTEXT));
			return (LRESULT) hBrushDeactive;
		}
		return 0;

  case WM_SETFOCUS :
    SetFocus(hEditIn);
    return 0;

	case WM_DESTROY :
    old_screen_off();
		DeleteObject(hBrushDeactive);        
		DeleteObject(hFontDeactive);        
		PostQuitMessage (0) ;
    //exit(0);
    return 0;

  case WM_CHAR:
    GetcBuf = wParam;
    SetEvent(KeyPress);
    if((wParam==3) && Processing==CLI_PROCESSING) 
      IShallDie = TRUE;
    break;

  case WM_COMMAND:
    id = LOWORD(wParam);
#ifdef DO_MENU
    if(id>=MIN_MENU && id<=MAX_MENU) gui_run_menu(id);
#endif
    return 0;
  }

  return DefWindowProc (hwnd, iMsg, wParam, lParam) ;
}

pict_info *get_pict_info(pict_info **list,HWND hwnd,bool remove){
  pict_info *curr = *list,*prev = NULL;
  while(curr && curr->hwnd != hwnd){
    prev = curr;
    curr = curr->next;
  }
  if(remove && curr){
    if(prev) prev->next = curr->next;
    else *list = curr->next;
  }
  return curr;
}

LRESULT CALLBACK PictProc (HWND hwnd, UINT iMsg, WPARAM wParam, LPARAM lParam){
  HDC hdc;
  HWND hpict;
  int step=0,low;

  if(iMsg==WM_VSCROLL || iMsg==WM_HSCROLL)
  {
    low = LOWORD(wParam);
    if(low==SB_THUMBPOSITION)
    {
      if(iMsg==WM_VSCROLL){
        yScroll = HIWORD(wParam);
        SetScrollPos(hwnd,SB_VERT,yScroll,TRUE);
      } else {
        xScroll = HIWORD(wParam);
        SetScrollPos(hwnd,SB_HORZ,xScroll,TRUE);
      }
    } else {
      if(low==SB_PAGEDOWN || low==SB_PAGEUP)
      {
        step = (int)(((iMsg==WM_VSCROLL)? MemHeight:MemWidth)/5);
        if(low==SB_PAGEUP) step = -step;
      } 
      else if(low==SB_LINEUP) step = -charheight;
      else if(low==SB_LINEDOWN) step = charheight;
      if(iMsg==WM_VSCROLL){
        yScroll = int_max(0,int_min(MemHeight-paneheight,yScroll+step));
        SetScrollPos(hwnd,SB_VERT,yScroll,TRUE);
      } else {
        xScroll = int_max(0,int_min(MemWidth-panewidth,xScroll+step));
        SetScrollPos(hwnd,SB_HORZ,xScroll,TRUE);
      }
    }
    RedrawPict(NULL);
    InvalidateRect(hPict,NULL,TRUE);
    UpdateWindow(hPict);
  }
  else switch(iMsg){
  case WM_SHOWSCROLL:
    if(lParam)
    {
      SetScrollPos(hwnd,SB_VERT,yScroll,FALSE);
      SetScrollPos(hwnd,SB_HORZ,xScroll,FALSE);
      SetScrollRange(hwnd,SB_VERT,0,int_max(MemHeight-1,paneheight-1),TRUE);
      SetScrollRange(hwnd,SB_HORZ,0,int_max(MemWidth-1,panewidth-1),TRUE);
    }
    ShowScrollBar(hwnd,SB_BOTH,lParam);
    return 0;

  case WM_LBUTTONDOWN:
    if(wParam & MK_SHIFT){
      hpict = CreateWindow (szPictName,         // window class name
						  "SPR Graphic",     // window caption
						  WS_OVERLAPPEDWINDOW,     // window style
	            CW_USEDEFAULT,           // initial x position
		          CW_USEDEFAULT,           // initial y position
			        493,                      // initial x size  -- could be 1024
				      512,                      // initial y size    -- could be 556
					    NULL,                    // parent window handle
						  NULL,                    // window menu handle
	            (HINSTANCE) GetWindowLong(hwnd,GWL_HINSTANCE),  // program instance handle
			        NULL) ;		             // creation parameters
      ShowWindow (hpict, 1) ;
      UpdateWindow (hpict) ;
    } else {
      hdc = GetDC(hwnd);
	    SetMapMode(hdc,MM_ISOTROPIC);
	    SetWindowExtEx(hdc,512,512,NULL);
	    SetViewportExtEx(hdc,clientwidth/2-2,-(int) (charheight*1.5-clientheight),NULL);
	    SetViewportOrgEx(hdc,0,0,NULL);
      MouseCoords.x = LOWORD(lParam);
      MouseCoords.y = HIWORD(lParam);
      MouseButton = 1;
      DPtoLP(hdc,&MouseCoords,1);
      MouseCoords.x = xScroll + MouseCoords.x;
      MouseCoords.y = MemHeight - yScroll - MouseCoords.y;
      if(Processing==CLI_COMMAND_OR_CLICK_INPUT){
        Processing = CLI_PROCESSING;
        SetEvent(ClickOrEnter);
      } else {
        SetEvent(MouseClick);
      }
	    ReleaseDC(hwnd,hdc);
    }
    return 0;

  case WM_RBUTTONDOWN:
    hdc = GetDC(hwnd);
	  SetMapMode(hdc,MM_ISOTROPIC);
	  SetWindowExtEx(hdc,512,512,NULL);
	  SetViewportExtEx(hdc,clientwidth/2-2,-(int) (charheight*1.5-clientheight),NULL);
	  SetViewportOrgEx(hdc,0,0,NULL);
    MouseCoords.x = LOWORD(lParam);
    MouseCoords.y = HIWORD(lParam);
    MouseButton = 3;
    DPtoLP(hdc,&MouseCoords,1);
    MouseCoords.x = xScroll + MouseCoords.x;
    MouseCoords.y = yScroll + 512 - MouseCoords.y;
    if(Processing==CLI_COMMAND_OR_CLICK_INPUT){
      Processing = CLI_PROCESSING;
      SetEvent(ClickOrEnter);
    } else {
      SetEvent(MouseClick);
    }
	  ReleaseDC(hwnd,hdc);
    return 0;
  }

  return DefWindowProc (hwnd, iMsg, wParam, lParam) ;
}
  
LRESULT CALLBACK PictWndProc (HWND hwnd, UINT iMsg, WPARAM wParam, LPARAM lParam){
  HDC hdc;
  PAINTSTRUCT ps ;
  pict_info *pinfo;

  switch(iMsg){
  case WM_CREATE:
    pinfo = (pict_info*) AM_MALLOC(pict_info);
    pinfo->next = pict_list;
    pict_list = pinfo;
    pinfo->width = 0;
    pinfo->height = 0;
    pinfo->hwnd = hwnd;

		hdc = GetDC(hwnd);
	  SetMapMode(hdc,MM_ISOTROPIC);
	  SetWindowExtEx(hdc,MemWidth,MemHeight,NULL);
	  SetViewportExtEx(hdc,pinfo->width,pinfo->height,NULL);
  	pinfo->hmemdc = CreateCompatibleDC(hdc);
		pinfo->hbitmap = CreateCompatibleBitmap(hdc,MemWidth,MemHeight);

	  SelectObject(pinfo->hmemdc,pinfo->hbitmap);
 	  BitBlt(pinfo->hmemdc,0,0,MemWidth,MemHeight,hPictDC,0,0,SRCCOPY);
    //SaveBitmap(hPictDC,hBitmap,"bmp.txt");

		ReleaseDC(hwnd,hdc);
    return 0;

  case WM_PAINT:
    pinfo = get_pict_info(&pict_list,hwnd,FALSE);
    hdc = BeginPaint (hwnd, &ps) ;
	  SetMapMode(hdc,MM_ISOTROPIC);
	  SetWindowExtEx(hdc,MemWidth,MemHeight,NULL);
	  SetViewportExtEx(hdc,pinfo->width,pinfo->height,NULL);
 	  BitBlt(hdc,0,0,MemWidth,MemHeight,pinfo->hmemdc,0,0,SRCCOPY);
    EndPaint (hwnd, &ps) ;
    return 0 ;

  case WM_SIZE:
    pinfo = get_pict_info(&pict_list,hwnd,FALSE);
    pinfo->width = LOWORD(lParam);
		pinfo->height = HIWORD(lParam);
    return 0;

  case WM_DESTROY:
    pinfo = get_pict_info(&pict_list,hwnd,TRUE);
  	DeleteDC(pinfo->hmemdc);
	  DeleteObject(pinfo->hbitmap);
    AM_FREE(pinfo,pict_info);
    return 0;
  }

  return DefWindowProc (hwnd, iMsg, wParam, lParam) ;
}

LRESULT CALLBACK EditProcProc (HWND hwnd, UINT iMsg, WPARAM wParam, LPARAM lParam){
  LRESULT ret;

  if(iMsg==WM_PAINT){
    ret=CallWindowProc(/*(int (__stdcall *)(void))*/ OldEditProcProc,hwnd,iMsg,wParam,lParam);
    SetEvent(SwitchedStatus);
    return ret;
  }
  else if(iMsg==WM_CHAR){
    GetcBuf = wParam;
    SetEvent(KeyPress);
    if((wParam==3) && Processing==CLI_PROCESSING) 
      IShallDie = TRUE;
  }
	return CallWindowProc(/*(int (__stdcall *)(void))*/ OldEditProcProc,hwnd,iMsg,wParam,lParam);
}

void refresh_command_history(){
  last_command = (last_command+1)%COMHISTSIZE;
  curr_command = (last_command+1)%COMHISTSIZE;
  strncpy(command_hist[last_command],EditInBuf,LINESIZE);
  strcpy(command_hist[curr_command],"");
}

LRESULT CALLBACK EditInProc (HWND hwnd, UINT iMsg, WPARAM wParam, LPARAM lParam)
{
  if(iMsg==WM_KEYDOWN){
    if(wParam==13 && Processing!=CLI_PROCESSING){
		  if(GetWindowText(hEditIn,(LPTSTR) EditInBuf,LINESIZE) 
          || Processing==CLI_SPEC_INPUT || Processing==CLI_KEY_INPUT){
        if(Processing==CLI_COMMAND_INPUT || Processing==CLI_COMMAND_OR_CLICK_INPUT){
          if(!strcmp(EditInBuf,"clear")){
            CaretPos = NextLine = NumOutLines = 0;
            strcpy(OutBuf,"");
            UpdateEditOut();
            SetWindowText(hEditIn,NULL);
          } else {
            refresh_command_history();
            if(Processing==CLI_COMMAND_INPUT)
              SetEvent(PressedEnter);
            else
              SetEvent(ClickOrEnter);
          }
        } else {
          SetEvent(PressedEnter);
        }
        Processing = CLI_PROCESSING;
  		  return 0;
		  }

    } else if(wParam==VK_UP){
      curr_command = curr_command-1;
      if(curr_command<0) curr_command+=COMHISTSIZE;
      SetWindowText(hEditIn,command_hist[curr_command]);
      return 0;

    } else if(wParam==VK_DOWN){
      curr_command = (curr_command+1)%COMHISTSIZE;
      if(curr_command==((last_command+2)%COMHISTSIZE)) 
        curr_command = (last_command+1)%COMHISTSIZE;
      SetWindowText(hEditIn,command_hist[curr_command]);
      return 0;
    }
	}
  else if(iMsg==WM_CHAR){
    GetcBuf = wParam;
    SetEvent(KeyPress);
    if(Processing==CLI_KEY_INPUT){
      sprintf(EditInBuf,"%c",EditInBuf[1]);
      SetWindowText(hEditIn,EditInBuf);
    }
    if((wParam==3) && Processing==CLI_PROCESSING) 
      IShallDie = TRUE;
  }

	return CallWindowProc(OldEditInProc,hwnd,iMsg,wParam,lParam);
}

void init_OFN(HWND hwnd){
  OFN.lStructSize = sizeof(OPENFILENAME);
  OFN.hwndOwner = hwnd;
  OFN.hInstance = NULL;
  OFN.lpstrCustomFilter = NULL;
  OFN.nMaxCustFilter = OFN.nFilterIndex = OFN.nFileOffset = OFN.nFileExtension = 0;
  OFN.lpstrFileTitle = (char*)malloc(sizeof(char)*50);
  OFN.nMaxFileTitle = 50;
  OFN.lpstrFile = (char*)malloc(sizeof(char)*100);
  OFN.nMaxFile = 100;
  OFN.lpstrInitialDir = NULL;
  OFN.lpstrTitle = NULL;
  OFN.lCustData = 0L;
  OFN.lpfnHook = NULL;
  OFN.lpTemplateName = NULL;
  /*These are good variables to set with the calling function*/
  OFN.Flags = 0;
  OFN.lpstrFilter = NULL;
  OFN.lpstrDefExt = NULL;
  strcpy(OFN.lpstrFileTitle,"");
  strcpy(OFN.lpstrFile,"");
}


void UpdateEditOut(){
	/*Change this to set text to only the portion of OutBuf that's in 
		the window if this is too slow...*/
	SetWindowText(hEditOut,OutBuf);
	SendMessage(hEditOut,EM_LINESCROLL,(WPARAM) 0,(LPARAM) NumOutLines-NumClientLines+2);
  FlushFlag = FALSE;
}

int fprintf(FILE *fil,const char *str,...){
  va_list args;
  char outstr[LINESIZE];
  char errbuf[100];
  int ret;

  va_start(args,str);
  if(fil!=stdout){
    ret=vfprintf(fil,str,args);
  } else {
	  ret=vsprintf(outstr,str,args);
  }
  va_end(args);

  if(ret<0){
    sprintf(errbuf,"fprintf() failed for unknown reasons.\nit may be that what you are printing exceeds the %d character buffer.",LINESIZE);
    my_error(errbuf);
  }
  if(fil==stdout) ret = screenprint(outstr);
  return ret;
}

int printf(const char *str,...){
	va_list args;
	char outstr[LINESIZE];
  char errbuf[100];
  int ret;

	va_start(args,str);
	ret=vsprintf(outstr,str,args);	
	va_end(args);
  if(ret<0){
    sprintf(errbuf,"printf() failed for unknown reasons.\nIt may be that what you are printing exceeds the %d character buffer.",LINESIZE);
    my_warning(errbuf);
  }
  return screenprint(outstr);
}

int screenprint(char *outstr){
  char temp[BUFSIZE];
	int i,nextline=0,len = strlen(outstr);
	
  if(LogFile) fprintf(LogFile,outstr);

	/*Change all '\n's to "\r\n"s */
	for(i=0;i<len;i++){
		if(outstr[i]=='\n' && (i==0 || outstr[i-1]!='\r')){
			len++;
			if(len>=LINESIZE) my_error("Buffer overflow in (f)printf!");
			strcpy(temp,outstr+i);
			strcpy(outstr+i+1,temp);
			outstr[i]='\r';
			if(!nextline) nextline = i+2;
			NumOutLines++;
		}
	}

	/*Delete previous lines until there's enough to append the new string*/
	while((CaretPos+len) >= BUFSIZE){
    if(NextLine){
		  memcpy(temp,OutBuf+NextLine,sizeof(char)*(CaretPos-NextLine+1));
		  memcpy(OutBuf,temp,sizeof(char)*(CaretPos-NextLine+1));
		  CaretPos -= NextLine;
		  NextLine = 0;
		  for(i=0;i<(CaretPos-1) && !NextLine;i++){
			  if(OutBuf[i]=='\r' && OutBuf[i+1]=='\n') NextLine = i+2;
		  }
    } else {
      len = 0;
      strcpy(outstr,"");
    }
	}

	/*Append the new string*/
	strcpy(OutBuf+CaretPos,outstr);
	CaretPos+=len;
	if(!NextLine) NextLine = nextline;

	if(nextline || FlushFlag) UpdateEditOut();

	return len;
}

void gui_flush(FILE *fil){
  if(fil==stdout) FlushFlag = TRUE;
}

char gui_getchar(void){
  char c;
  EnterCriticalSection(&critsect);
  ResetEvent(KeyPress);
  WaitForSingleObject(KeyPress,INFINITE);
  c = GetcBuf;
  LeaveCriticalSection(&critsect);
  return ((c=='\r')? '\n':c);
}

char gui_getc(FILE *s)
{
  if ( s == stdin )
    return gui_getchar();
  else
    return fgetc(s);
}

char* gui_getstring(FILE *s){
  char *str;
  ResetEvent(PressedEnter);
  if(s==stdin){
    WaitForSingleObject(PressedEnter,INFINITE);
    str = mk_copy_string(EditInBuf);
  } else {
    str = AM_MALLOC_ARRAY(char,LINESIZE);
    fgets(str,LINESIZE,s);
  }
  return str;
}

void gui_prompt(char *str,int wait_type,bool fullsize){
	XSplit = (fullsize)? clientwidth : (charwidth*((str)? strlen(str):0));
  SetWindowText(hEditIn,NULL);
  SetWindowText(hEditProc,str);
  MoveWindow(hEditProc,0,(int) (clientheight-charheight*1.5),XSplit,(int) (charheight*1.5),TRUE);
  MoveWindow(hEditIn,XSplit,(int) (clientheight-charheight*1.5),clientwidth,(int) (charheight*1.5),TRUE);
  SetFocus(hEditIn);
  if(str) ShowWindow(hEditProc,SW_SHOWNORMAL);
  ShowWindow(hEditIn,SW_SHOWNORMAL);
  Processing = wait_type;
}

void gui_unprompt(void){
  SetFocus(hPict);
  Processing = CLI_PROCESSING;
  ShowWindow(hEditIn,SW_HIDE);
  SetWindowText(hEditProc,"Processing....");
  MoveWindow(hEditProc,0,(int) (clientheight-charheight*1.5),clientwidth,(int) (charheight*1.5),TRUE);
	MoveWindow(hEditIn,0,(int) (clientheight-charheight*1.5),clientwidth,(int) (charheight*1.5),TRUE);
}

void gui_quit(void){
  int i;
  for(i=0;i<main_argc;i++) free(main_argv[i]);
  free(main_argv);
  free(OFN.lpstrFileTitle);
  free(OFN.lpstrFile);
  if(LogFile) fclose(LogFile);
  SetWindowText(hEditProc,"Quit.");
  DeleteCriticalSection(&critsect);
//  PostQuitMessage(0);
//  _endthread();
  exit(0);
}

void gui_link_command_to_menu(char *name,int id){
#ifdef DO_MENU
  int i;

  i = id-MIN_MENU;
  if(!(MenuCommand[i].id)){
    strcpy(MenuCommand[i].name,name);
    MenuCommand[i].id = id;
  }
#endif
}

void gui_run_menu(int id){
#ifdef DO_MENU
  menu_item *m;
  char *str = NULL;

  if(Processing==CLI_COMMAND_OR_CLICK_INPUT){
    m = &MenuCommand[id-MIN_MENU];
    if(m->id && m->name && (int)strlen(m->name)>0)
    {
      strcpy(EditInBuf,m->name);
      refresh_command_history();
      ResetEvent(FinishedCommand);
      Processing = CLI_PROCESSING;
      SetEvent(ClickOrEnter);
    }
  }
#endif
}

#else /*now if !PC_MVIS_PLATFORM*/

#ifndef UNIX_PLATFORM

int fprintf(FILE *fil,const char *str,...){
  va_list args;
  int ret;

  va_start(args,str);
  if(fil!=stdout){
    ret=vfprintf(fil,str,args);
  } else {
	  ret=vprintf(str,args);
    if(LogFile) vfprintf(LogFile,str,args);
  }
  va_end(args);

  if(ret<0){
    char errbuf[128];
    sprintf(errbuf,"fprintf() failed for unknown reasons.\nIt may be that what you are printing exceeds the %d character buffer.",LINESIZE);
    my_error(errbuf);
  }
  return ret;
}

int printf(const char *str,...){
	va_list args;
  int ret;

	va_start(args,str);
	ret=vprintf(str,args);
  if(LogFile) vfprintf(LogFile,str,args);
	va_end(args);

  if(ret<0){
    char errbuf[100];
    sprintf(errbuf,"printf() failed for unknown reasons.");
    my_warning(errbuf);
  }
  return ret;
}
#endif /* ! UNIX_PLATFORM */

void gui_flush(FILE *fil){
  fflush(fil);
}

char gui_getchar(void){
  return getchar();
}

char gui_getc(FILE *s){
  return fgetc(s);
}

char* gui_getstring(FILE *s){
  char *str = AM_MALLOC_ARRAY(char,LINESIZE);
  fgets(str,LINESIZE,s);
  return str;
}

void gui_prompt(char *prompt,int wait_type,bool fullsize){
  printf("%s",prompt);
}

void gui_unprompt(void){
}

void gui_link_command_to_menu(char *name,int id){
}

void gui_quit(void){
  if(LogFile) fclose(LogFile);
  exit(0);
}

#endif  /*PC_MVIS_PLATFORM*/

void gui_log(char *filename){
  if(LogFile) fclose(LogFile);
  if(filename){
    LogFile = fopen(filename,"a");
    if(!LogFile){
      char errbuf[100];
      sprintf(errbuf,"The file %s could not be opened for logging",filename);
      my_warning(errbuf);
    }
  } else LogFile = NULL;
}

