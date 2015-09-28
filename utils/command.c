/*
   File:        command.c
   Author:      Andrew W. Moore
   Created:     Mon Jun  2 19:58:44 EDT 1997
   Description: Utilities on datsets

   Copyright 1998, SPR
*/

#include "command.h"
#include "adgui.h"

#ifdef PC_MVIS_PLATFORM
extern int Processing;
extern HANDLE PressedEnter,SwitchedStatus,ClickOrEnter,FinishedCommand;
extern char EditInBuf[];
extern POINT MouseCoords;
extern int MouseButton;
extern OPENFILENAME OFN;
#endif

command *mk_command_from_string_array(string_array *sa)
{
  command *c = AM_MALLOC(command);

  c -> line = mk_copy_string_array(sa);
  c -> args = mk_copy_string_array(sa);
  string_array_remove(c->args,0);
  c -> s  = (string_array_size(sa) > 0) ? mk_copy_string(string_array_ref(sa,0)) : NULL;
  c -> a1 = (string_array_size(sa) > 1) ? mk_copy_string(string_array_ref(sa,1)) : NULL;
  c -> a2 = (string_array_size(sa) > 2) ? mk_copy_string(string_array_ref(sa,2)) : NULL;
  c -> a3 = (string_array_size(sa) > 3) ? mk_copy_string(string_array_ref(sa,3)) : NULL;
  c -> a4 = (string_array_size(sa) > 4) ? mk_copy_string(string_array_ref(sa,4)) : NULL;
  c -> attnums = NULL;

  return c;
}

void rebuild_command_from_string_array(command *c,string_array *sa){
  command *new_c = mk_command_from_string_array(sa);

  free_string_array(c->line);
  free_string_array(c->args);
  if(c->s) free_string(c->s);
  if(c->a1) free_string(c->a1);
  if(c->a2) free_string(c->a2);
  if(c->a3) free_string(c->a3);
  if(c->a4) free_string(c->a4);
  if(c->attnums) free_ivec(c->attnums);

  c->line = mk_copy_string_array(new_c->line);
  c->args = mk_copy_string_array(new_c->args);
  c->s = (new_c->s)? mk_copy_string(new_c->s):NULL;
  c->a1 = (new_c->a1)? mk_copy_string(new_c->a1):NULL;
  c->a2 = (new_c->a2)? mk_copy_string(new_c->a2):NULL;
  c->a3 = (new_c->a3)? mk_copy_string(new_c->a3):NULL;
  c->a4 = (new_c->a4)? mk_copy_string(new_c->a4):NULL;
  c->attnums = NULL;

  free_command(new_c);
}

command *mk_command_from_string(char *s)
{
  string_array *sa = mk_broken_string(s);
  command *c = mk_command_from_string_array(sa);
  free_string_array(sa);
  return c;
}

#ifdef PC_MVIS_PLATFORM

command *mk_command_from_line(char *prompt)
{
  command *c;
  string_array *sa = NULL;
  Processing = CLI_COMMAND_OR_CLICK_INPUT;

  while(Processing!=CLI_PROCESSING){
    MouseButton = 0;
    gui_prompt(NULL,CLI_COMMAND_OR_CLICK_INPUT,FALSE);
    WaitForSingleObject(ClickOrEnter,INFINITE);
    if(MouseButton)
      sprintf(EditInBuf,"%s-click %d %d",(MouseButton==1)? "left":"right"
                                            ,MouseCoords.x,MouseCoords.y);
    //else printf("\nYou entered: %s\n",EditInBuf);
    sa = mk_broken_string(EditInBuf);
    if ( sa == NULL ){
      Processing = CLI_COMMAND_OR_CLICK_INPUT;
    }
    else if ( string_array_size(sa) < 1 ){
      free_string_array(sa);
      Processing = CLI_COMMAND_OR_CLICK_INPUT;
    }
  }
  ResetEvent(SwitchedStatus);
  gui_unprompt();

  WaitForSingleObject(SwitchedStatus,INFINITE);
  if( (string_array_size(sa)>1) && !strcmp(string_array_ref(sa,1),"=") ) // this added by Artur to allow for 'x = y' version of  'define x y'
  {
	  string_array_set(sa,1,string_array_ref(sa,0));
	  string_array_set(sa,0,"define");
  }
  c = mk_command_from_string_array(sa);
  free_string_array(sa);
  return c;
}

#else

command *mk_command_from_line(char *prompt)
{
  command *c;
  string_array *sa = NULL;
  while ( sa == NULL )
  {
    printf("%s",prompt);
    sa = mk_string_array_from_line(stdin);
    if ( sa == NULL )
    {
      /* skip */
    }
    else if ( string_array_size(sa) < 1 )
    {
      free_string_array(sa);
      sa = NULL;
    }
  }
  c = mk_command_from_string_array(sa);
  free_string_array(sa);
  return c;
}

#endif

char *command_arg_ref(command *c,int arg)
{
  char *result = NULL;
  if ( arg < string_array_size(c->line) )
    result = string_array_ref(c->line,arg);
  return result;
}

int command_num_tokens(command *c)
{
	return (c->line)? string_array_size(c->line):0;
}

char *basic_string_from_command(command *c,char *key)
{
  char *result = NULL;
  int i;

  if (c)
  {
    for ( i = 1 ; result == NULL && i < command_num_tokens(c)-1 ; i++ )
      if ( arg_is(c,i,key) ) result = command_arg_ref(c,i+1);
  }

  return result;
}

char *mk_possibly_empty_string_from_user(char *prompt)
{
  char buff[1000];
  int i = 0;
  int c = ' ';
  string_array *sa;
  char *string;

#ifdef PC_MVIS_PLATFORM
  gui_prompt(prompt,CLI_SPEC_INPUT,FALSE);
  string = gui_getstring(stdin);
  strcpy(buff,string);
  free_string(string);
  gui_unprompt();
#else
  printf("%s> ",prompt);
  while ( i < 999 && c != EOF && c != '\n' )
  {
    c = getc(stdin);
    if ( c != EOF && c != '\n' )
    {
      buff[i] = (char) c;
      i += 1;
    }
  }
  buff[i] = '\0';
#endif

  /* Following lines remove whitespace from start and end */
  sa = mk_broken_string(buff);
  string = mk_string_from_string_array(sa);
  free_string_array(sa);
  return string;
}

bool is_a_boolean(char *s)
{
  bool result;
  if ( s == NULL )
    result = FALSE;
  else
  {
    char c = s[0];
    if ( c >= 'a' && c <= 'z' ) c += 'A' - 'a';
    result = (c == '0' || c == '1' || c == 'N' || c == 'Y' || c == 'T' || c == 'F');
  }
  return result;
}

bool atob(char *s)
{
  char c;
  if ( !is_a_boolean(s) )
    my_error("atob passed non-boolean string");
  c = s[0];
  if ( c >= 'a' && c <= 'z' ) c += 'A' - 'a';
  return c == '1' || c == 'Y' || c == 'T';
}

bool validate_arg(char *arg,int argtype,double argmin,double argmax,string_array *argparams){
  bool ret = (arg!=NULL);
  int iarg;
  double darg;
  if(argtype==ARG_FILE_SAVE) argtype = ARG_FILE_OPEN;
  else if(argtype==ARG_STRING_SET) argtype = ARG_STRING;
  else if(argtype>=ARG_INT_LIST) argtype -= ARG_INT_LIST;

  if(arg) switch(argtype){
  case ARG_FILE_OPEN:
    ret = TRUE;
    break;
  case ARG_STRING:
    ret = !argparams || find_index_in_string_array(argparams,arg)!=-1;
    break;
  case ARG_INT:
    ret = is_a_number(arg);
    if(ret) {
      iarg = atoi(arg);
      ret = (argmin==0 && argmax==0) || (iarg>=(int)argmin && iarg<=(int)argmax);
    }
    break;
  case ARG_REAL:
    ret = is_a_number(arg);
    if(ret) {
      darg = atof(arg);
      ret = (argmin==0 && argmax==0) || (darg>=argmin && darg<=argmax);
    }
    break;
  case ARG_BOOL:
    ret = is_a_boolean(arg);
    break;
  default:
    my_error("Cannot validate: invalid argument type");
  }
  return ret;
}

/* Message may be NULL, which means it's ignored.

   If non-null, then if we end up prompting the user,
   we prompt them with the message instead of the key */
char *mk_string_from_command_with_message(command *c,char *key,char *message,
                                          char *default_string)
{
  char *x = basic_string_from_command(c,key);
  char *result = NULL;
  char prompt[1000];

  if ( x != NULL )
    result = mk_copy_string(x);
  else
  {
    sprintf(prompt,"%s [default %s]",(message==NULL)?key:message,default_string);
    result = mk_possibly_empty_string_from_user(prompt);    

    if ( result == NULL || result[0] == '\0' )
    {
      if ( result != NULL ) free_string(result);
      result = mk_copy_string(default_string);
    }
  }

  return result;
}

char *mk_string_from_command(command *c,char *key,char *default_string)
{
  return mk_string_from_command_with_message(c,key,NULL,default_string);
}

int int_from_command(command *c,char *key,int default_val,int lo,int hi)
{
  bool ok = FALSE;
  char *base = basic_string_from_command(c,key);
  int result = -7777;
  char *string = (base==NULL) ? NULL : mk_copy_string(base);
  char *prompt = NULL;

  default_val = int_min(hi,default_val);
  default_val = int_max(lo,default_val);

  while ( !ok )
  {
    ok = validate_arg(string,ARG_INT,lo,hi,NULL);
    if(ok)
		  result = atoi(string);
    else
    {
      if ( string != NULL )
        printf("The value of %s should be an integer\n"
               "between %d and %d. %s is unsuitable.\n",key,lo,hi,string);
      
      if ( string != NULL ) free_string(string);
      prompt = mk_printf("Enter %s (%d to %d) [default %d]",key,lo,hi,default_val);
      string = mk_possibly_empty_string_from_user(prompt);
      if ( string[0] == '\0' )
      {
        ok = TRUE;
        result = default_val;
      }
    }
  }

  if ( string != NULL ) free_string(string);
  if (prompt != NULL) free_string(prompt);

  return result;
}

double double_from_command(command *c,char *key,double default_val,
			   double lo,double hi)
{
  bool ok = FALSE;
  char *base = basic_string_from_command(c,key);
  double result = -7777.7;
  char *string = (base==NULL) ? NULL : mk_copy_string(base);
  char prompt[100];

  default_val = real_min(hi,default_val);
  default_val = real_max(lo,default_val);

  while ( !ok )
  {
    ok = validate_arg(string,ARG_REAL,lo,hi,NULL);
    if ( ok )
      result = atof(string);
    else
    {
      if ( string != NULL )
        printf("The value of %s should be a real number\n"
               "between %g and %g. %s is unsuitable.\n",key,lo,hi,string);
      
      if ( string != NULL ) free_string(string);
      sprintf(prompt,"Enter %s (%g to %g) [default %g] ",key,lo,hi,default_val);
      string = mk_possibly_empty_string_from_user(prompt);
      if ( string[0] == '\0' )
      {
        ok = TRUE;
        result = default_val;
      }
    }
  }

  if ( string != NULL ) free_string(string);

  return result;
}

bool bool_from_command(command *c,char *key,bool default_val)
{
  bool ok = FALSE;
  char *base = basic_string_from_command(c,key);
  bool result = FALSE;
  char *string = (base==NULL) ? NULL : mk_copy_string(base);
  char prompt[100];

  while ( !ok )
  {
    ok = validate_arg(string,ARG_BOOL,0,0,NULL);
    if(ok) 
      result = atob(string);
    else
    {
      if ( string != NULL )
        printf("The value of %s should be a truth-value.\n"
               "%s is unsuitable.\n",key,string);
      
      if ( string != NULL ) free_string(string);
      sprintf(prompt,"Enter choice %s [default %s] ",key,(default_val)?"y":"n");
      string = mk_possibly_empty_string_from_user(prompt);
      if ( string[0] == '\0' )
      {
        ok = TRUE;
        result = default_val;
      }
    }
  }

  if ( string != NULL ) free_string(string);

  return result;
}


bool is_ok_number_between(command *c,int arg,double lo,double hi)
{
  char *name = command_arg_ref(c,arg);
  bool result = FALSE;
  if ( name != NULL && is_a_number(name) )
  {
    double x = atof(name);
    result = x >= lo && x <= hi;
  }
  if ( !result )
  {
    printf("I needed item %d on the line to be a number between\n"
           "%g and %g inclusive. ",arg,lo,hi);
    if ( name == NULL )
      printf("But the item was missing\n");
    else
      printf("But %s is not acceptable.\n",name);
  }
  return result;
}

int filter_length(string_array *sa){
  int i,len=1;
  for(i=1;i<string_array_size(sa);i++)
    len += ((int)strlen(string_array_ref(sa,i))+1);
  return len;
}

void free_filter_from_string_array(char *filter,string_array *sa){
  AM_FREE_ARRAY(filter,char,filter_length(sa));
}

char *mk_filter_from_string_array(string_array *sa){
  int i,len = filter_length(sa);
  char *s = AM_MALLOC_ARRAY(char,len);
  s[len-1] = '\0';
  len = 0;

  for(i=1;i<string_array_size(sa);i++){
    strcpy(s+len,string_array_ref(sa,i));
    len += ((int)strlen(string_array_ref(sa,i))+1);
  }
  return s;
}

char *mk_string_from_default_arg(int type,double def,string_array *params){
  char *ret = NULL;
  if(type==ARG_STRING_SET) type=ARG_STRING;
  else if(type>=ARG_INT_LIST) type-=ARG_INT_LIST;
  if(type==ARG_STRING || type==ARG_FILE_OPEN || type==ARG_FILE_SAVE)
    ret = mk_copy_string(params? string_array_ref(params,(int)def):"");
  else {
    char *buff = NULL;
    switch(type){
    case ARG_INT:
      buff=mk_printf("%d",(int)def); break;
    case ARG_REAL:
      buff=mk_printf("%f",def); break;
    case ARG_BOOL:
      buff=mk_printf("%s",params? string_array_ref(params,(int)def):(def?"Yes":"No")); break;
    default:
      my_error("Trying to make default argument of invalid type");
    }
    ret = mk_copy_string(buff);
    free_string(buff);
  }
  return ret;
}

void *mk_arg_from_string(char *arg,int type){
  void *ret;
  string_array *sa;
  int i;

  if(type==ARG_STRING || type==ARG_FILE_OPEN || type==ARG_FILE_SAVE)
    ret = mk_copy_string(arg);
  else if(type==ARG_STRING_LIST || type==ARG_STRING_SET)
    ret = mk_broken_string(arg);
  else switch(type){
    case ARG_INT:
      ret = AM_MALLOC(int);
      *(int*)ret = atoi(arg);
      break;
    case ARG_REAL:
      ret = AM_MALLOC(double);
      *(double*)ret = atof(arg);
      break;
    case ARG_BOOL:
      ret = AM_MALLOC(bool);
      *(bool*)ret = atob(arg);
      break;
    case ARG_INT_LIST:
      sa = mk_broken_string(arg);
      ret = mk_ivec(0);
      for(i=0;i<string_array_size(sa);i++) add_to_ivec((ivec*)ret,atoi(string_array_ref(sa,i)));
      free_string_array(sa);
      break;
    case ARG_REAL_LIST:
      sa = mk_broken_string(arg);
      ret = mk_dyv(0);
      for(i=0;i<string_array_size(sa);i++) add_to_dyv((dyv*)ret,atof(string_array_ref(sa,i)));
      free_string_array(sa);
      break;
    default:
      my_error("Trying to build arg of invalid type\n");
  }
  return ret;
}

#ifdef PC_MVIS_PLATFORM
typedef struct arg_dlg{
  string_array *values;
  int type;
  double min,max,def;
  string_array *params;
  char name[50];
  char caption[100];
  char text[200];
  bool hit_cancel;
} arg_dlg;

void *mk_arg_dlg(int type,char *name,double min,double max,double def,string_array *params){
  arg_dlg *arg = AM_MALLOC(arg_dlg);
  arg->values = NULL;
  arg->type = (type==ARG_BOOL)? ARG_STRING:type;
  arg->min = min;
  arg->max = max;
  arg->def = (max!=min)? min(max,max(def,min)):def;
  arg->params = params;
  arg->hit_cancel = FALSE;
  if(type==ARG_STRING_SET) sprintf(arg->name,"strings");
  else switch(type%ARG_INT_LIST){
  case ARG_INT:
    sprintf(arg->name,"an integer"); break;
  case ARG_REAL:
    sprintf(arg->name,"a real"); break;
  case ARG_STRING:
    sprintf(arg->name,"a string"); break;
  case ARG_BOOL:
    sprintf(arg->name,"a boolean"); 
    break;
  }
  sprintf(arg->caption,"Choose %s%s",(type==ARG_STRING_SET)?"a set of ":"",name? name:arg->name);
  if(max!=min){
    if(type==ARG_INT)
      sprintf(arg->text,"Please choose %s value between %d and %d%s%s.",arg->name,
                        (int)min,(int)max,name? " to represent ":"",name? name:"");
    else 
      sprintf(arg->text,"Please choose %s value between %f and %f%s%s.",arg->name,
                        min,max,name? " to represent ":"",name? name:"");
  } else {
    if(type==ARG_BOOL && name && 
        (count_occurences_of_char_in_string(name,' ')>5 || strchr(name,'?')))
      sprintf(arg->text,"%s",name);
    else
      sprintf(arg->text,"Please choose %s value%s%s.",arg->name,
                        name? " to represent ":"",name? name:"");
  }
  return arg;
}

void free_arg_dlg(arg_dlg *arg){
  if(arg->values) free_string_array(arg->values);
  AM_FREE(arg,arg_dlg);
}

BOOL CALLBACK process_arg(HWND hDlg,UINT iMsg,WPARAM wParam,LPARAM lParam){
  static char buff[500];
  static double value;
  static int spin_val;
  static HWND hw,hwspin;
  arg_dlg *arg = (arg_dlg*) DlgInfo;
  int i;
  char *temp;

  switch(iMsg){
  case WM_INITDIALOG:
    value = arg->def;
    temp = mk_string_from_default_arg(arg->type,arg->def,arg->params);
    strcpy(buff,temp);
    free_string(temp);
    SetWindowText(hDlg,arg->caption);
    SetWindowText(GetDlgItem(hDlg,IDC_ARG_DESCRIPTION),arg->text);

    switch((arg->type)%ARG_INT_LIST){
    case ARG_INT:
      hw = GetDlgItem(hDlg,IDC_EDIT_INT);
      hwspin = GetDlgItem(hDlg,IDC_SPIN_INT);
      ShowWindow(hwspin,SW_SHOWNORMAL); 
      SetWindowText(hw,buff);
      break;
    case ARG_REAL:
      hw = GetDlgItem(hDlg,IDC_EDIT_ARG);
      SetWindowText(hw,buff);
      break;
    case ARG_STRING:
      if(arg->params){
        hw = GetDlgItem(hDlg,IDC_DROP_STRING);
        for(i=0;i<string_array_size(arg->params);i++) 
          SendMessage(hw,CB_ADDSTRING,0,(LPARAM)string_array_ref(arg->params,i));
        SendMessage(hw,CB_SETCURSEL,(WPARAM)(int)value, 0L);
      } else {
        hw = GetDlgItem(hDlg,IDC_EDIT_ARG);
      }
      break;
    }

    if(arg->type>=ARG_INT_LIST){
      ShowWindow(GetDlgItem(hDlg,IDNEXT),SW_SHOWNORMAL);
      ShowWindow(GetDlgItem(hDlg,IDFINISHED),SW_SHOWNORMAL);
      ShowWindow(GetDlgItem(hDlg,IDCANCEL2),SW_SHOWNORMAL);
      ShowWindow(GetDlgItem(hDlg,IDOK),SW_HIDE);
      ShowWindow(GetDlgItem(hDlg,IDCANCEL),SW_HIDE);
      EnableWindow(GetDlgItem(hDlg,IDOK),FALSE);
    } else {
      EnableWindow(GetDlgItem(hDlg,IDNEXT),FALSE);
    }
    ShowWindow(hw,SW_SHOWNORMAL); 
    SetFocus(hw);
    return TRUE;
    break;

  case WM_VSCROLL:
    if(LOWORD(wParam)==4){
      GetWindowText(hw,(LPTSTR)buff,50);
      value = atoi(buff) + ((HIWORD(wParam)<spin_val || HIWORD(wParam)==0)? 1:-1);
      if(arg->min!=arg->max) value = min(arg->max,max(arg->min,value));
      sprintf(buff,"%d",(int)value);
      SetWindowText(hw,buff);
      spin_val = HIWORD(wParam);
    }
    break;

  case WM_COMMAND:
    if(LOWORD(wParam)==IDOK || LOWORD(wParam)==IDNEXT || LOWORD(wParam)==IDFINISHED){
      GetWindowText(hw,(LPTSTR)buff,50);
      if(!validate_arg(buff,arg->type,arg->min,arg->max,arg->params)) {
        char warning[100];
        sprintf(warning,"\'%s\' is not a valid value",buff);
        my_warning(warning);
      } else {
        if(!arg->values) arg->values = mk_string_array(0);
        add_to_string_array(arg->values,buff);
        if(LOWORD(wParam)==IDNEXT){
          gui_dialog(IDD_ARG,process_arg);
          if(arg->hit_cancel){
            arg->hit_cancel = FALSE;
            string_array_remove_last_element(arg->values);
            SetFocus(hDlg);
            return FALSE;
          }
        }
        EndDialog(hDlg,0);
        return TRUE;
      }
    } else if(LOWORD(wParam)==IDCANCEL || LOWORD(wParam)==IDCANCEL2){
      arg->hit_cancel = TRUE;
      EndDialog(hDlg,0);
      return TRUE;
    }
  }
  return FALSE;
}

BOOL CALLBACK process_argset(HWND hDlg,UINT iMsg,WPARAM wParam,LPARAM lParam){
  arg_dlg *arg = (arg_dlg*) DlgInfo;
  HWND hlist = GetDlgItem(hDlg,IDC_LIST);
  int i,max;
  char buff[1000];
  switch(iMsg){
  case WM_INITDIALOG:
    SetWindowText(hDlg,arg->caption);
    for(i=0;i<string_array_size(arg->params);i++)
      SendMessage(hlist,LB_ADDSTRING,0,(LPARAM)string_array_ref(arg->params,i));
    return TRUE;
  case WM_COMMAND:
    switch(LOWORD(wParam)){
    case IDCANCEL:
      arg->hit_cancel = TRUE;
      EndDialog(hDlg,0);
      return TRUE;
    case IDOK:
      max = SendMessage(hlist,LB_GETCOUNT,0,0);
      if(!arg->values) arg->values = mk_string_array(0);
      for(i=0;i<max;i++){
        if(SendMessage(hlist,LB_GETSEL,i,0)){
          SendMessage(hlist,LB_GETTEXT,i,(LPARAM)buff);
          add_to_string_array(arg->values,buff);
        }
      }
      EndDialog(hDlg,1);
      return TRUE;
    }
  }
  return FALSE;
}

void* mk_arg_string_from_invalid_command(command *c,int argnum,char *argname,int argtype,
                                         double argmin,double argmax,double argdefault,
                                         string_array *argparams,bool validation){
  char *arg = NULL;

  gui_prompt("",CLI_PROCESSING,TRUE);

  if(argtype==ARG_FILE_OPEN || argtype==ARG_FILE_SAVE)
  {
    char *ext = mk_copy_string(string_array_ref(argparams,0));
    char *filter = mk_filter_from_string_array(argparams);
    char *title = mk_copy_string((argname)? argname:((argtype==ARG_FILE_OPEN)? "Open File":"Save File"));
    OFN.Flags = OFN_HIDEREADONLY | ((argtype==ARG_FILE_SAVE)? OFN_OVERWRITEPROMPT : 0);
    OFN.lpstrDefExt = ext;
    OFN.lpstrTitle = title;
    OFN.lpstrFilter = filter;
    if((argtype==ARG_FILE_OPEN)? GetOpenFileName(&OFN) : GetSaveFileName(&OFN)){
      gui_prompt((argtype==ARG_FILE_OPEN)? "Loading....":"Saving....",CLI_PROCESSING,TRUE);
      arg = mk_copy_string(OFN.lpstrFile);
    }
    free_string(title);
    free_string(ext);
    free_filter_from_string_array(filter,argparams);
  } 
  else 
  {
    DlgInfo = mk_arg_dlg(argtype,argname,argmin,argmax,argdefault,argparams);
    if(argtype==ARG_BOOL){
      bool barg = MessageBox(hWnd,((arg_dlg*)DlgInfo)->text,((arg_dlg*)DlgInfo)->caption
                                 ,MB_YESNOCANCEL|((argdefault)? MB_DEFBUTTON1:MB_DEFBUTTON2));
      if(barg!=IDCANCEL) arg = mk_copy_string((barg==IDYES)? "Yes":"No");
    } else {
      if(argtype==ARG_STRING_SET) gui_dialog(IDD_ARGSET,process_argset);
      else gui_dialog(IDD_ARG,process_arg);
      if(((arg_dlg*)DlgInfo)->values){
        arg = mk_string_from_string_array(((arg_dlg*)DlgInfo)->values);
      }
    }
    free_arg_dlg((arg_dlg*)DlgInfo);
    gui_prompt("Processing....",CLI_PROCESSING,TRUE);
  }
  return arg;
}

#else

void* mk_arg_string_from_invalid_command(command *c,int argnum,char *argname,int argtype,
                                         double argmin,double argmax,double argdefault,
                                         string_array *argparams,bool validation){
  char *arg = NULL;
  string_array *sa = NULL;
  char *buff;
  char *my_argname = mk_copy_string(argname? argname:"value");
  bool ok = FALSE;

  if(argtype>=ARG_INT_LIST) sa = mk_string_array(0);

  while (!ok)
  {
    if (arg && arg[0]=='\0')
    {
      free_string(arg);
      arg = mk_string_from_default_arg(argtype,argdefault,argparams);
    }
    if(arg) ok = !validation || validate_arg(arg,argtype,argmin,argmax,argparams);

    if (!ok || (sa && strcmp(arg,"done")))
    {
      if ( arg ) {
        if(!ok) {
          printf("%s is an invalid value.\n",arg);
        } else {
          add_to_string_array(sa,arg);
        }
        free_string(arg);
      }
      if (argtype==ARG_STRING || argtype==ARG_STRING_LIST || argtype==ARG_STRING_SET || argtype==ARG_BOOL){
        buff=mk_printf("Enter %s ",my_argname);
        if(argparams) 
	{
	  char *new_buff = mk_printf("%s [default %s]",buff,string_array_ref(argparams,(int)argdefault));
	  free_string(buff);
	  buff = new_buff;
	}
        if(argtype>=ARG_INT_LIST) 
        {
	  char *new_buff = mk_printf("%s, or 'done' to finish",buff);
	  free_string(buff);
	  buff = new_buff;
	}
      }
      else
        buff=mk_printf("Enter %s (%g to %g%s) [default %g] ",my_argname,argmin,argmax,
                        (argtype>=ARG_INT_LIST)?", or 'done' to finish":"",argdefault);
      arg = mk_possibly_empty_string_from_user(buff);
      free_string(buff);
      ok = FALSE;
    }
  }
  if(sa)
  {
    free_string(arg);
    arg = mk_string_from_string_array(sa);
    free_string_array(sa);
  }
  free_string(my_argname);

  return arg;
}
#endif

char *mk_arg_string_from_valid_command(command *c,int argnum,char *argname,int argtype,
                                       double argmin,double argmax,double argdef,
                                       string_array *argparams,bool validation){
  char *arg = NULL;
  string_array *sa = NULL;
  int num;
  bool ok = FALSE;

  if(argname)
  {
    num = find_index_in_string_array(c->args,argname)+2;
    if(num!=1 && (!validation || validate_arg(command_arg_ref(c,num),argtype,argmin,argmax,argparams))){
      argnum = num;
      ok = TRUE;
    }
  }
  if(!ok && (!validation || validate_arg(command_arg_ref(c,argnum),argtype,argmin,argmax,argparams)))
    ok = TRUE;
  if(!ok) my_error("mk_arg_string_from_valid_command was called on an invalid command");

  arg = command_arg_ref(c,argnum);
  if(argtype>=ARG_INT_LIST) 
  {
    sa = mk_string_array(0);
    while(ok && strcmp(arg,"done"))
    {
      add_to_string_array(sa,arg);
      arg = command_arg_ref(c,++argnum);
      ok = validate_arg(arg,argtype,argmin,argmax,argparams);
    }
  }

  if(sa){
    arg = mk_string_from_string_array(sa);
    free_string_array(sa);
  } else arg = mk_copy_string(arg);

  return arg;    
}

void* mk_arg_from_command(command *c,int argnum,char *argname,int argtype,
                          double argmin,double argmax,double argdefault,
                          string_array *params,bool validation){
  void *ret = NULL;
  string_array *argparams = NULL;
  string_array *sa,*sa_new;
  char *arg = NULL;
  int i;
  bool insert_arg = TRUE;

  if(params) argparams = mk_copy_string_array(params);
  else if(argtype==ARG_BOOL) argparams = mk_string_array_2("No","Yes");
  if((argtype==ARG_INT || argtype==ARG_INT_LIST || argtype==ARG_REAL || argtype==ARG_REAL_LIST) 
      && (argmin!=0 || argmax!=0)){
    if(argdefault>argmax) argdefault = argmax;
    else if(argdefault<argmin) argdefault = argmin;
  }
      
  if(argname){
    arg = basic_string_from_command(c,argname);
    if(arg)
      if(validation && !validate_arg(arg,argtype,argmin,argmax,argparams)) arg = NULL;
  }
  if(!arg && argnum!=-1){
    arg = command_arg_ref(c,argnum);
    if(arg)
      if(validation && !validate_arg(arg,argtype,argmin,argmax,argparams)) arg = NULL;
  }

  if(arg) {
    insert_arg = FALSE;
    arg = mk_arg_string_from_valid_command(c,argnum,argname,argtype,argmin,argmax,argdefault,argparams,validation);
  } else {
    arg = (char*)mk_arg_string_from_invalid_command(c,argnum,argname,argtype,argmin,argmax,argdefault,argparams,validation);
  }

  if(arg){
    if(insert_arg){
      sa = mk_copy_string_array(c->line);
      sa_new = mk_broken_string(arg);
      if(argnum==-1) {
        argnum=1+find_index_in_string_array(sa,argname);
        if(!argnum) {
          add_to_string_array(sa,"done");
          add_to_string_array(sa,argname);
          argnum = string_array_size(sa);
        }
      }
      for(i=string_array_size(sa_new)-1;i>=0;i--)
        insert_in_string_array(sa,argnum,string_array_ref(sa_new,i));
      rebuild_command_from_string_array(c,sa);
      free_string_array(sa);
      free_string_array(sa_new);
    }
    ret = mk_arg_from_string(arg,argtype);
    free_string(arg);
  }

  if(argparams) free_string_array(argparams);
  return ret;
}

void free_arg(void *arg,int argtype){
  if(arg) {
    if(argtype==ARG_STRING || argtype==ARG_FILE_OPEN || argtype==ARG_FILE_SAVE)
      free_string((char*)arg);
    else if(argtype==ARG_STRING_LIST || argtype==ARG_STRING_SET)
      free_string_array((string_array*)arg);
    else switch(argtype){
      case ARG_INT:
        AM_FREE(arg,int); break;
      case ARG_REAL:
        AM_FREE(arg,double); break;
      case ARG_BOOL:
        AM_FREE(arg,bool); break;
      case ARG_INT_LIST:
        free_ivec((ivec*)arg); break;
      case ARG_REAL_LIST:
        free_dyv((dyv*)arg); break;
      default:
        my_error("Trying to free invalid argument type");
    }
  }
}

bool arg_exists_ok(command *c,int arg)
{
  char *name = command_arg_ref(c,arg);
  bool result = name != NULL;
  if ( !result )
    printf("You are missing at least one argument on the line\n");
  return result;
}

bool arg_is(command *c,int argnum,char *name)
{
  char *arg = command_arg_ref(c,argnum);
  bool result;

  if ( name == NULL )
    result = arg == NULL;
  else if ( arg == NULL )
    result = FALSE;
  else
    result = caseless_eq_string(name,arg);

  return result;
}

void free_command(command *c)
{
  free_string_array(c->line);
  free_string_array(c->args);
  if ( c -> s != NULL ) free_string(c->s);
  if ( c -> a1 != NULL ) free_string(c->a1);
  if ( c -> a2 != NULL ) free_string(c->a2);
  if ( c -> a3 != NULL ) free_string(c->a3);
  if ( c -> a4 != NULL ) free_string(c->a4);
  if ( c -> attnums != NULL ) free_ivec(c->attnums);
  AM_FREE(c,command);
}

char *Help_string = "Undefined";
int Menu_ID = 0;
bool Acted = FALSE;
void *GlobalEnv = NULL;

typedef struct helper
{
  string_array *names;
  string_array *briefs;
  string_array *details;
} helper;

helper *mk_empty_helper()
{
  helper *he = AM_MALLOC(helper);
  /* Chasing a memory leak? It's not your fault! All you need to do is
     put free_global_helper() at the end of your main.c before your
     final am_malloc_report()*/
  he -> names = mk_string_array(0);
  he -> briefs = mk_string_array(0);
  he -> details = mk_string_array(0);
  return he;
}

/* Finds hindex such that name == helper_name(he,hindex).
   Returns -1 if no such name */
int helper_name_to_hindex(helper *he,char *name)
{
  return index_in_sosarray(he->names,name);
}

bool help_exists(helper *he,char *name)
{
  int hindex = helper_name_to_hindex(he,name);
  return hindex >= 0;
}

void add_helper_string(helper *he,char *name,char *brief,char *detail)
{
  if ( !help_exists(he,name) )
  {
    int index = find_sosarray_insert_index(he->names,name);
    insert_in_string_array(he->names,index,name);
    insert_in_string_array(he->briefs,index,brief);
    insert_in_string_array(he->details,index,detail);
  }
}

int find_index_of_char_in_string(char *s,char c)
{
  int result = -1;
  int i;
  for ( i = 0 ; result < 0 && s[i] != '\0' ; i++ )
  {
    if ( c == s[i] )
      result = i;
  }
  return result;
}

int find_index_in_string(char *s,char *stops)
{
  int result = -1;
  int i;
  for ( i = 0 ; result < 0 && s[i] != '\0' ; i++ )
  {
    if ( 0 <= find_index_of_char_in_string(stops,s[i]) )
      result = i;
  }
  return result;
}

/* Makes a new null-terminated string using character
   s[start] s[start+1] ... s[end-1]
   strlen(result) = end - start

   PRE: 0 <= start < end <= strlen(s) 
*/
char *mk_substring(char *s,int start,int end)
{
  char *result;
  int i;
  if ( start < 0 || start >= end || end > (int)strlen(s) )
    my_error("mk_substring bad start or end");

  result = AM_MALLOC_ARRAY(char,end-start+1); /* +1 for '\0' */

  for ( i = 0 ; i < end-start ; i++ )
    result[i] = s[start+i];
  result[end-start] = '\0';
  return result;
}

/* Just like mk_substring except returns an empty string ("")
   if start >= end, and makes start at least 0 and end at most
   strlen(s) before going into action. Always succeeds. */
char *mk_friendly_substring(char *s,int start,int end)
{
  char *result;

  start = int_max(0,start);
  end = int_min((int)strlen(s),end);

  if ( end <= start )
    result = mk_copy_string("");
  else
    result = mk_substring(s,start,end);

  return result;
}

char *mk_remove_chars(char *string,char *seppers)
{
  int temp_size = (int)strlen(string) + 1;
  char *temp = AM_MALLOC_ARRAY(char,temp_size);
  char *result = NULL;
  int string_index;
  int temp_index = 0;
  for ( string_index = 0 ; string[string_index] != '\0' ; string_index++ )
  {
    char c = string[string_index];
    if ( find_index_of_char_in_string(seppers,c) < 0 )
    {
      temp[temp_index] = c;
      temp_index += 1;
    }
  }
  temp[temp_index] = '\0';
  result = mk_copy_string(temp);
  AM_FREE_ARRAY(temp,char,temp_size);
  return result;
}

int helper_size(helper *he)
{
  return string_array_size(he->names);
}

char *helper_name_ref(helper *he,int hindex)
{
  return string_array_ref(he->names,hindex);
}

char *helper_brief_ref(helper *he,int hindex)
{
  return string_array_ref(he->briefs,hindex);
}

char *helper_detail_ref(helper *he,int hindex)
{
  return string_array_ref(he->details,hindex);
}

bool help_string_matches(helper *he,int hindex,char *key)
{
  bool result = caseless_eq_string(helper_name_ref(he,hindex),key);
  return result;
}

int find_index_or_end_in_string(char *string,char *stops)
{
  int index = find_index_in_string(string,stops);
  if ( index < 0 )
    index = (int) strlen(string);
  return index;
}

void parse_helper_string(char *help_string,char **r_name,char **r_brief,char **r_detail)
{
  int space_index = find_index_or_end_in_string(help_string," \n");
  int return_index = find_index_or_end_in_string(help_string,"\n");
  char *first_token = mk_friendly_substring(help_string,0,space_index);

  *r_name = mk_remove_chars(first_token,":");
  *r_brief = mk_friendly_substring(help_string,space_index+1,return_index);
  *r_detail = mk_friendly_substring(help_string,return_index,(int)strlen(help_string));

  free_string(first_token);
}

void parse_and_add_helper_string(helper *he,char *help_string)
{
  char *name;
  char *brief;
  char *detail;

  parse_helper_string(help_string,&name,&brief,&detail);
  add_helper_string(he,name,brief,detail);
  
  free_string(name);
  free_string(brief);
  free_string(detail);
}

void free_helper(helper *he)
{
  free_string_array(he->names);
  free_string_array(he->briefs);
  free_string_array(he->details);
  AM_FREE(he,helper);
}

void display_help_string(char *help_string,bool one_line)
{
  if ( !one_line )
  {
    printf("%s",help_string);
    printf("\n");
  }
  else
  {
    int i = 0;
    while ( help_string[i] != '\n' && help_string[i] != '\0' )
    {
      printf("%c",help_string[i]);
      i += 1;
    }
    printf("\n");
  }
}

helper *Global_helper = NULL;

helper *get_global_helper()
{
  if ( Global_helper == NULL )
  {
    Global_helper = mk_empty_helper();
  }
  return Global_helper;
}

void free_global_helper()
{
  if ( Global_helper != NULL )
  {
    free_helper(Global_helper);
    Global_helper = NULL;
  }
}

void command_malloc_report()
{
  if ( Global_helper != NULL )
  {
    fprintf(stdout,"#\n#IMPORTANT IMPORTANT IMPORTANT ************\n"
	           "# The cached data structures used by command.c\n"
	           "# are still allocated.\n"
	           "# Call free_global_helper() [prototype in command.h] if\n"
                   "# you wish to deallocate them.\n");
  }
}

#define LONG_MODE    0
#define ONELINE_MODE 1
#define NAME_MODE    2

void display_one_helper(FILE *s,helper *he,int hindex,int display_mode)
{
  fprintf(s,"%s",helper_name_ref(he,hindex));
  if ( display_mode == NAME_MODE )
    fprintf(s,"\n");
  else
  {
    fprintf(s,": %s\n",helper_brief_ref(he,hindex));
    if ( display_mode == LONG_MODE )
      fprintf(s,"%s\n",helper_detail_ref(he,hindex));
  }
}

string_array *mk_helper_names(helper *he,ivec *hindexes)
{
  int size = ivec_size(hindexes);
  string_array *names = mk_string_array(size);
  int i;
  for ( i = 0 ; i <size ; i++ )
  {
    int hindex = ivec_ref(hindexes,i);
    char *name = helper_name_ref(he,hindex);
    string_array_set(names,i,name);
  }
  return names;
}

void display_helper_names(FILE *s,helper *he,ivec *hindexes)
{
  string_array *sa = mk_helper_names(he,hindexes);
  display_names(s,sa);
  free_string_array(sa);
}

void display_helper(FILE *s,helper *he,ivec *hindexes,int display_mode)
{
  if ( display_mode == NAME_MODE )
    display_helper_names(s,he,hindexes);
  else
  {
    int i;
    for ( i = 0 ; i < ivec_size(hindexes) ; i++ )
    {
      int hindex = ivec_ref(hindexes,i);
      display_one_helper(s,he,hindex,display_mode);
    }
  }
}

void save_helper_prolog(FILE *s)
{
  fprintf(s,"<HEAD>\n");
  fprintf(s,"<TITLE>Schenley Park Research: Available Commands</TITLE>\n");
  fprintf(s,"</HEAD>\n");
  fprintf(s,"<BODY>\n");
  fprintf(s,"<H1>Schenley Park Research: Available Commands</H1>\n");
  fprintf(s,"<P>\n");
}

void save_helper_epilog(FILE *s)
{
  fprintf(s,"<P>\n");
  fprintf(s,"</BODY>\n");
}

void html_string_out(FILE *s,char *string)
{
  int i;
  for ( i = 0 ; string[i] != '\0' ; i++ )
  {
    char c = string[i];
    if ( c == '<' )
      fprintf(s,"&#60;");
    else if ( c == '>' )
      fprintf(s,"&#62;");
    else
      fprintf(s,"%c",c);
  }
}

void save_one_helper_to_html(FILE *s,helper *he,int hindex)
{
  char *name = helper_name_ref(he,hindex);
  fprintf(s,"<hr>\n");
  fprintf(s,"<H2><A NAME=\"%s\">%s</A></H2>\n",name,name);
  fprintf(s,"<P>\n");
  fprintf(s,"<H4>%s : ",name);
  html_string_out(s,helper_brief_ref(he,hindex));
  fprintf(s,"</H4>\n");
  fprintf(s,"<PRE>\n");
  html_string_out(s,helper_detail_ref(he,hindex));
  fprintf(s,"\n</PRE>\n");
}

void save_helper_to_html(FILE *s,helper *he)
{
  int i;
  save_helper_prolog(s);
  for ( i = 0 ; i < helper_size(he) ; i++ )
  {
    char *name = helper_name_ref(he,i);
    fprintf(s,"<A HREF=#%s>%s</A>%s\n",
           name,name,(i==helper_size(he)-1)?"":", ");
  }

  fprintf(s,"<hr>\n");

  for ( i = 0 ; i < helper_size(he) ; i++ )
    save_one_helper_to_html(s,he,i);
  save_helper_epilog(s);
}

bool hindex_matches_key(helper *he,int hindex,char *key)
{
  bool result;
  if ( caseless_eq_string(key,"all") )
    result = TRUE;
  else 
  {
    char *name = helper_name_ref(he,hindex);
    result = string_pattern_matches(key, name);
    /*
    if ( key[0] == name[0] && key[1] == '\0' )
      result = TRUE;
    else
      result = caseless_eq_string(key,name);
      */
  }
  return result;
}

ivec *mk_hindexes_matching_key(helper *he,char *key)
{
  ivec *hindexes = mk_ivec(0);
  int i;
  for ( i = 0 ; i < helper_size(he) ; i++ )
  {
    if ( hindex_matches_key(he,i,key) )
      add_to_ivec(hindexes,i);
  }
  return hindexes;
}

void give_help(FILE *s,helper *he,char *key)
{
  ivec *hindexes = mk_hindexes_matching_key(he,key);
  int size = ivec_size(hindexes);
  if ( size == 0 )
    fprintf(s,"I can't find any help matching \"%s\"\n",key);
  else
  {
    int display_mode = (size == 1) ? LONG_MODE : ( size < 15 ) ? ONELINE_MODE : NAME_MODE;
    display_helper(s,he,hindexes,display_mode);
  }
  free_ivec(hindexes);
}

void display_global_helper(char *key)
{
  helper *he = get_global_helper();
  give_help(stdout,he,key);
}

void add_to_global_helper(char *help_string)
{
  helper *he = get_global_helper();
  parse_and_add_helper_string(he,help_string);
}

void generic_dispatch(void *env,command *c,char *name,
                      void (*generic_dispatch)(void *env,command *c,bool act))
{
  bool give_help = caseless_eq_string(c->s,"help");
  bool build_menu = eq_string(c->s,"rebuild_menu");

  if(build_menu)
  {
    Menu_ID = 0;
    generic_dispatch(env,c,FALSE);
    if(Menu_ID) gui_link_command_to_menu(name,Menu_ID);
    Acted = TRUE;
  }
  else if ( give_help )
  {
    generic_dispatch(env,c,FALSE);
    add_to_global_helper(Help_string);
    Acted = TRUE;
  }
  else if ( caseless_eq_string(c->s,name) )
  {
    generic_dispatch(env,c,TRUE);
    Acted = TRUE;
  }
}

command *Next_Command = NULL;
bool Stop_Commands = FALSE;
int Num_Iterations = 0;
typedef struct dispatch_node{
  command *c;
  bool stop;
  int num_its;
} dispatch_node;

dispatch_node *mk_dispatch_node(){
  dispatch_node *dn = AM_MALLOC(dispatch_node);
  dn->c = Next_Command;
  dn->stop = Stop_Commands;
  dn->num_its = Num_Iterations;
  return dn;
}

void free_dispatch_node(void *vdn)
{
  dispatch_node *dn = (dispatch_node *)vdn;
  free_command(dn->c);
  AM_FREE(dn,dispatch_node);
}

void multiple_dispatch(void *env,command *c,char *name,
                       void (*local_dispatch)(void *env,command *c,bool act),
                       void dispatch_command(void *env,command *c))
{
  bool act = caseless_eq_string(c->s,name);
  static generic_array *ga = NULL;
  if(act){
    Next_Command = NULL;
    Stop_Commands = FALSE;
    Num_Iterations = 0;
    if(!ga) ga = mk_empty_generic_array(free_dispatch_node,NULL,NULL);
  }
  generic_dispatch(env,c,name,local_dispatch);
  if(act){
    while(Next_Command){
      dispatch_node *dn;
      add_pointer_to_generic_array(ga,mk_dispatch_node());
      dispatch_command(env,Next_Command);
      dn = (dispatch_node *)generic_array_ref(ga,generic_array_size(ga)-1);
      Stop_Commands = dn->stop;
      Num_Iterations = dn->num_its+1;
      Next_Command = NULL;
      generic_array_remove(ga,generic_array_size(ga)-1);
      if(!Stop_Commands) local_dispatch(env,c,TRUE);
    }
    if(generic_array_size(ga)==0){
      free_generic_array(ga);
      ga = NULL;
    }
  }
}

void dispatch_quit(void *env,command *c,bool act)
{
  Help_string = "quit: end program\n"
                "Warning---please make sure you've saved any results you\n"
                "need before you type quit\n";
#ifdef DO_MENU
  Menu_ID = ID_FILE_QUIT;
#endif
  if ( act )
    printf("Bye bye!\n");
}

bool Comment_on = FALSE;
char Session_name[1000];
int Graphic_number = 0;

bool comment_is_on()
{
  return Comment_on;
}

char *comment_session_name()
{
  if ( !comment_is_on() )
    my_error("comment_session_name(): only valid if command is on");
  return Session_name;
}

void switch_comment_on(char *session_name)
{
  FILE *s = stdout;
  Comment_on = TRUE;
  sprintf(Session_name,"%s",session_name);
  Graphic_number = 1;
  fprintf(s,"\n");
  fprintf(s,"\n");
  fprintf(s,"\n");
  fprintf(s,"\\documentstyle[psfig,fullpage]{article}\n");
  fprintf(s,"\n");
  fprintf(s,"\\newcommand{\\bigfig}[1]{\\mbox{\\begin{minipage}{0.6\\textwidth}\n");
  fprintf(s,"\\centerline{\\psfig{file=#1,width=\\textwidth}}\n");
  fprintf(s,"\\end{minipage}}\\ \\\\\n");
  fprintf(s,"\\noindent}\n");
  fprintf(s,"\n");
  fprintf(s,"\\title{\\bf Example session using an SPR application}\n");
  fprintf(s," \n");
  fprintf(s,"\\author{{\\bf My Name Here}\\\\ \n");
  fprintf(s,"My Company\\\\ \n");
  fprintf(s,"My Street\\\\\n");
  fprintf(s,"My City\\\\\n");
  fprintf(s,"My Phone\\\\\n");
  fprintf(s,"My Email}\n");
  fprintf(s,"\n");
  fprintf(s,"\\begin{document} \n");
  fprintf(s,"\n");
  fprintf(s,"\\maketitle \n");
  fprintf(s,"\n");
  fprintf(s,"\\begin{verbatim}\n");
}

void switch_comment_off()
{
  FILE *s = stdout;
  Comment_on = FALSE;
  fprintf(s,"  \n");
  fprintf(s,"\n");
  fprintf(s,"\\end{verbatim}\n");
  fprintf(s,"\n");
  fprintf(s,"\\end{document}\n");
  fprintf(s,"\n");
  fprintf(s,"\n");
}

char *mk_line_from_user(char *prompt)
{
  char *string;
  gui_prompt(prompt,CLI_SPEC_INPUT,FALSE);
  string = gui_getstring(stdin);
  gui_unprompt();
  return string;
}

string_array *mk_paragraph_from_user()
{
  bool finished = FALSE;
  string_array *lines = mk_string_array(0);
  while ( !finished )
  {
    char *prompt = mk_printf("Line %d of your paragraph (. to end)> ",
                             1 + string_array_size(lines));
    char *line = mk_line_from_user(prompt);
    if ( eq_string(line,".") )
      finished = TRUE;
    else
      add_to_string_array(lines,line);

    free_string(line);
    free_string(prompt);
  }
  return lines;
}

void output_comment(FILE *s,string_array *lines)
{
  int i;
  if ( comment_is_on() )
    fprintf(s,"\\end{verbatim}\n\n\\begin{quote}\n");
  for ( i = 0 ; i < string_array_size(lines) ; i++ )
    fprintf(s,"%s\n",string_array_ref(lines,i));
  if ( comment_is_on() )
    fprintf(s,"\\end{quote}\n\n\\begin{verbatim}\n");
}

void dispatch_comment(void *env,command *c,bool act)
{
  Help_string = "comment <words> : Add a comment to a demo session.\n"
   "This command is almost certainly useless to you unless you\n"
   "are planning on saving the text of your session to a file, and wish\n"
   "to add a comment to the reader of the session.\n"
   "\n"
   "comment on <sessionname>: This will switch on some internal flags that will\n"
   "            add special annotations to all output so that it can easily\n"
   "            be turned into a tutorial. Graphics will be saved to files\n"
   "            with names like andrew12.ps if andrew was given as the\n"
   "            session name.\n"
   "\n"
   "comment off: Switch off the above.\n"
   "\n"
   "comment <anything else>: Will allow the user to type in commentary about\n"
   "                         what's going on to be used by the comment command.\n"
   "         Important: TO STOP ENTERING COMMENT, PUT A SINGLE . ON A LINE.\n";
    
  if ( act ) 
  {
    bool on = arg_is(c,1,"on");
    bool off = arg_is(c,1,"off");
    if ( on )
    {
      if ( arg_exists_ok(c,2) )
      {
        char *session_name = c->a2;
        switch_comment_on(session_name);
      }
    }
    else if ( off )
      switch_comment_off();
    else
    {
      char *firstline = mk_string_from_string_array(c->args);
      string_array *lines = mk_paragraph_from_user();
      insert_in_string_array(lines,0,firstline);
      output_comment(stdout,lines);
      free_string_array(lines);
      free_string(firstline);
    }
  }
}

void dispatch_pause(void *env,command *c,bool act){
  if(act)
  {
    wait_for_key();
  }
}

void dispatch_log(void *env,command *c,bool act){
  Help_string = "log <filename/off> : Enable or disable logging of all text output to a file\n"
                "This command takes all text output and appends it into a file.\n"
                "(Note that appending allows the user to add to a log file at any time.)\n"
                "Typing \'log off\' at any time will disable logging.";
  if(act)
  {
    char *arg = c->a1;
    string_array *params = mk_string_array_5("txt","Text Files (*.TXT)","*.txt"
                                                  ,"All Files (*.*)","*.*");
    arg = (char *)mk_arg_from_command(c,1,"File for logging",ARG_FILE_OPEN,0,0,0,params,1);
    if(arg)
    {
      if(caseless_eq_string(arg,"off")) gui_log(NULL);
      else gui_log(arg);
    }
    free_string_array(params);
    free_arg(arg,ARG_FILE_OPEN);
  }
}

void dispatch_batch(void *env,command *c,bool act)
{
  Help_string = "batch <filename>: execute the series of commands in <filename>\n"
                "The commands are executed in just the same way as if you\n"
                "were to type them into the command line manually here. Note\n"
                "that all output is sent to this display as usual. Note too\n"
                "that the execution of the commands continues blindly even\n"
                "if there are errors along the way.\n\n"
                "The batch file should simply have one command on each line.\n"
                "You can leave blank lines if you like. Lines beginning\n"
               "with a # character will be treated as comments and ignored.\n";

#ifdef DO_MENU
  Menu_ID = ID_FILE_BATCH;
#endif
  if ( act ) {
    void generic_batch(void *env,FILE *f);
    static FILE* f = NULL;
    if(Num_Iterations==0) {
      char *filename = NULL;
      string_array *sa = mk_string_array_5("BAT","Batch Files (*.BAT)","*.bat"
                                                ,"Text Files (*.TXT)","*.txt");
      add_to_string_array(sa,"All Files (*.*)");
      add_to_string_array(sa,"*.*");
      filename = (char*)mk_arg_from_command(c,1,"Batch file to execute",ARG_FILE_OPEN,0,0,0,sa,1);
      free_string_array(sa);
      if(filename){
        f = fopen(filename,"r");
        if(!f){
          char *buff=mk_printf("Can't open file %s to run batch commands\n",filename);
          my_warning(buff);
          free_string(buff);
        } else printf("Executing batch file %s:\n",filename);
      }
      else printf("You must supply a filename\n");
      free_arg(filename,ARG_FILE_OPEN);
    }
    if(f){
      generic_batch(env,f);
      if(Stop_Commands) fclose(f);
    }
    else Stop_Commands = TRUE;
    Acted = TRUE;
  }
}

void *Parent_Env = NULL;
void (*Parent_GC)(void *env,command *c) = NULL;

/* Returns TRUE if acted */
void basic_command(void *env,command *c)
{
  if ( arg_is(c,0,"help") )
    printf("(Help options: type help, or help all, or help <command>)\n");
  generic_dispatch(env,c,"quit",dispatch_quit);
  generic_dispatch(env,c,"comment",dispatch_comment);
  generic_dispatch(env,c,"log",dispatch_log);
  generic_dispatch(env,c,"pause",dispatch_pause);
  multiple_dispatch(Parent_Env,c,"batch",dispatch_batch,Parent_GC);
}

bool Already_saved_html = FALSE;

void execute_generic_command(void *env,command *c,
                             void (*generic_command)(void *env,command *c))
{
  Acted = FALSE;
  Parent_Env = env;
  Parent_GC = generic_command;

  generic_command(env,c);
  if ( !Acted )
    printf("*** Ignoring unrecognized command %s\n",c->s);

  if ( eq_string(c->s,"help") )
  {
    if ( !Already_saved_html )
    {
      char *helpfile = "help.html";
      FILE *s = fopen(helpfile,"w");
      if ( s == NULL )
        printf("I would have saved all help information to %s,\n"
               "except I don't have write permission\n",helpfile);
      else
      {
        void save_helper_to_html(FILE *s,helper *he);
        save_helper_to_html(s,Global_helper);
        fclose(s);
        printf("*** I have saved all the help information in %s\n",helpfile);
      }
      Already_saved_html = TRUE;
    }
    display_global_helper((c->a1==NULL)?"all":c->a1);
  }
}

void generic_batch(void *env,FILE *s)
{
  static bool oldsm;
  static int linenum;
  string_array *sa;
  if(Num_Iterations==0){
    oldsm = SuppressMessages;
    linenum = 0;
  }
  gui_prompt("Executing batch file....",CLI_PROCESSING,TRUE);
  SuppressMessages = TRUE;

  sa = mk_next_interesting_line(s,&linenum);

  if ( sa == NULL ){
    Stop_Commands = TRUE;
    SuppressMessages = oldsm;
    printf("\n----------------------------------------------------------\n"
           "Finished executing from batch file\n");
  }
  else
  {
		char *string; 
		if( (string_array_size(sa)>1) && !strcmp(string_array_ref(sa,1),"=") ) // this added by Artur to allow for 'x = y' version of  'define x y'
		{
			string_array_set(sa,1,string_array_ref(sa,0));
			string_array_set(sa,0,"define");
		}
		Next_Command = mk_command_from_string_array(sa);
		string = mk_string_from_string_array(sa);
    printf("\n----------------------------------------------------------\n"
           "BATCH line %d> %s\n",linenum,string);
		free_string(string);
    free_string_array(sa);
  }
}

void rebuild_menu(void *env,void (*generic_command)(void *env,command *c)){
  command *c;
  string_array *sa = mk_string_array_1("rebuild_menu");
  c = mk_command_from_string_array(sa);
  execute_generic_command(env,c,generic_command);
  free_command(c);
  free_string_array(sa);
}

void generic_cli(void *env,char *prompt,
                 void (*generic_command)(void *env,command *c),
		 int argc,char *argv[])
{
  char *batchfile = string_from_args("batch",argc,argv,NULL);

  bool finished = FALSE;
#ifdef DO_MENU
  rebuild_menu(env,generic_command);
#endif

  if ( batchfile != NULL )
  {
    if ( !file_exists(batchfile) )
      my_errorf("Can't run batchfile %s...no such file\n",batchfile);
    else
    {
      char *cs = mk_printf("batch %s",batchfile);
      command *c = mk_command_from_string(cs);
      printf("*** Before entering the interactive command line, I'll execute\n"
	     "*** from file %s mentioned in the program command line\n",
	     batchfile);

      execute_generic_command(env,c,generic_command);
      free_command(c);
      free_string(cs);
    }
  }

  while ( !finished )
  {
    command *c = mk_command_from_line(prompt);
    char *string = mk_string_from_string_array(c->line);

    if ( comment_is_on() )
    {
      if ( !eq_string(c->s,"comment") )
      {
        printf("\n----------------------------------------------------------\n");
        printf("\\end{verbatim}\n");
        printf("{\\large \\bf\n");
        printf("\\begin{verbatim}\n");
        printf("SPR> %s\n",string);
        printf("\\end{verbatim}}\n");
        printf("\\begin{verbatim}\n");
      }
    }
    else if(!strstr(string,"-click "))
      printf("\n----------------------------------------------------------\n"
             "SPR> %s\n",string);

    if ( caseless_eq_string(c->s,"quit") )
      finished = TRUE;
    else
      execute_generic_command(env,c,generic_command);

    free_string(string);
    free_command(c);
  }
}

void basic_cli(int argc,char *argv[])
{
  void *env = NULL;
  generic_cli(env,"basic> ",basic_command,argc,argv);
}

void basic_cli_main(int argc,char *argv[])
{
  basic_cli(argc,argv);
}

/* Useful if you are parsing a command and want to look for a certain key
   and then remove it. Returns TRUE if it found the key. */
bool find_and_remove(string_array *sa,char *key)
{
  int index = find_index_in_string_array(sa,key);
  bool result = index >= 0;
  if ( result ) 
    string_array_remove(sa,index);
  return result;
}

/* string_has_suffix("plop.ps",".ps") == TRUE
   string_has_suffix("plop.ps","ps") == TRUE
   string_has_suffix("plops",".ps") == FALSE
   string_has_suffix("plops","ps") == TRUE
   string_has_suffix("plop.ps.hello",".ps") == FALSE */
bool string_has_suffix(char *string,char *suffix)
{
  int suffix_length = (int) strlen(suffix);
  int string_length = (int) strlen(string);
  int matches = string_length >= suffix_length;
  int i;
  for ( i = 0 ; matches && i < suffix_length ; i++ )
    matches = string[string_length - suffix_length + i] == suffix[i];
  return matches;
}

/* Finds least index such that
     string_has_suffix(sa[index],suffix)
   and returns -1 if no such string exists in sa */
int find_index_of_string_with_suffix(string_array *sa,char *suffix)
{
  int result = -1;
  int i;
  for ( i = 0 ; result < 0 && i < string_array_size(sa) ; i++ )
  {
    if ( string_has_suffix(string_array_ref(sa,i),suffix) )
      result = i;
  }
  return result;
}

/* If there's a string in sa which has the given suffix, returns a
   copy of the leftmost (lowest indexed such string) and removes that
   one string from sa.

   If there's no such string, returns NULL and leaves sa unchanged */
char *mk_find_string_with_suffix_and_remove(string_array *sa,char *suffix)
{
  int index = find_index_of_string_with_suffix(sa,suffix);
  char *result = (index < 0) ? NULL : 
                 mk_copy_string(string_array_ref(sa,index));
  if ( index >= 0 )
    string_array_remove(sa,index);
  return result;
}

/* Plop */

