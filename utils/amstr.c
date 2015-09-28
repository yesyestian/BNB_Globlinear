
/*
   File:         amstr.c
   Author:       Andrew Moore
   Created:      June 25 1997 from amma.c by Frank Dellaert
   Description:  String related functions.

   Copyright (C) Andrew W. Moore, 1992
*/

#include "amstr.h"
#include <ctype.h>

char *make_copy_string(const char *s)
{
  char *newv;


  if ( s == NULL )
  {
    newv = NULL;
  }
  else
  {
    newv = AM_MALLOC_ARRAY(char,1+strlen(s));
    sprintf(newv,"%s",s);
  }

  return(newv);
}

char *mk_printf(char *fmt, ...)
{
  va_list ap;
  char buff[10000];

  va_start(ap, fmt);
  vsprintf(buff, fmt, ap);
  va_end(ap);
  return (mk_copy_string(buff));
}

/* This function returns a string which is forced to have
   the given extension.  If the string already has the
   extension, the string is copied and returned.  If it
   does not have the extension, then that extension is
   added without any other modification to the string
   or other extensions it may already have.  The extension
   should be given without the ".", and the dot will be
   added if the extension is added.
*/
char *mk_string_extension(char *s, char *ext)
{
  char *sub = s;
  char *dot_ext = am_malloc(strlen(ext)+2);
  char *res;
  bool has_ext = FALSE;

  dot_ext[0] = '.';
  strcpy(dot_ext+1,ext);
  while(sub)
  {
    sub = strstr(sub,dot_ext);
    if (sub && (strlen(sub) == strlen(dot_ext))) 
    {
      has_ext = TRUE;
      break;
    }
    if (sub) sub++;  /* make sure we keep progressing through the string */
  }
  if (has_ext) res = mk_copy_string(s);
  else
  {
    res = am_malloc(strlen(s) + strlen(dot_ext) + 1);
    strcpy(res,s);
    strcpy(res+strlen(s),dot_ext);
  }
  am_free_string(dot_ext);
  return res;
}

char *mk_copy_string(char *s)
{
  return(make_copy_string(s));
}

/* s is given in a possibly quoted format.  all delimiting quotes are removed.
   all backslash specified characters are replaced with their appropriate
   single characters.  There is one additional special case.  The string
   NULL is converted to the NULL pointer.  If the actual string NULL is wanted
   it should be specified as 'NULL' or "NULL"
   Modified 1-31-96 JS: single quotes are no longer valid delimiters
 */
char *mk_copy_quoted_string(char *s)
{
  int i,index;
  char *newv;
  char *resultv;

  if ((!s) || (!strcmp(s,"NULL"))) return NULL;

  newv = AM_MALLOC_ARRAY(char,1+strlen(s));

  for (i=0,index=0;i< (int) strlen(s);i++)
  {
    if (s[i] == '\"') continue;
    else if (s[i] == '\\')
    {
      if ((s[i+1] == '\"')||(s[i+1] == '\\')) 
      {
	newv[index++] = s[i+1]; 
	i++;
      }
      /* technically there shouldn't be an else, but we'll be graceful */
      else newv[index++] = s[i];
    }
    else newv[index++] = s[i];
  }
  newv[index] = '\0';
  resultv = mk_copy_string(newv);
  am_free(newv,1+strlen(s));
  return resultv;
}

/* inverts the mk_copy_quoted_string operation */
char *mk_copy_to_quoted_string(char *s)
{
  int i,size,index;
  char *newv;
  bool whitespace = FALSE;

  if (!s) return mk_copy_string("NULL");

  for (i=0,size=2;i< (int) strlen(s);i++)
  {
    if (isspace(s[i])) whitespace = TRUE;
    if ((s[i]=='\\')||(s[i]=='\"')) size+=2;
    else size++;
  }

  if (!whitespace) return mk_copy_string(s);

  newv = AM_MALLOC_ARRAY(char,1+size);
  index = 0;
  newv[index++] = '\"';
  for (i=0;i< (int) strlen(s);i++)
  {
    if ((s[i]=='\\')||(s[i]=='\"')) newv[index++] = '\\';
    newv[index++] = s[i];
  }
  newv[index++] = '\"';
  newv[index] = '\0';
  return newv;
}

char *mk_input_string(char *message)
{
  char buff[200];
  printf("%s > ",message);
  input_string("",buff,200);
  return(mk_copy_string(buff));
}

char *mk_downcase_string(char *s)
{
  char *result = mk_copy_string(s);
  int i;
  for ( i = 0 ; result[i] != '\0' ; i++ )
    if ( result[i] >= 'A' && result[i] <= 'Z' )
      result[i] = (char)(result[i] - 'A' + 'a');
  return(result);
}

char *mk_downcase_string_with_length(char *s,int n)
{
  char *result = mk_copy_string(s);
  int i;
  for ( i = 0 ; i < n ; i++ )
    if ( result[i] >= 'A' && result[i] <= 'Z' )
      result[i] = (char)(result[i] - 'A' + 'a');
  return(result);
}

char *mk_upcase_string(char *s)
{
  char *result = mk_copy_string(s);
  int i;
  for ( i = 0 ; result[i] != '\0' ; i++ )
    if ( result[i] >= 'a' && result[i] <= 'z' )
      result[i] = (char)(result[i] - 'a' + 'A');
  return(result);
}

char *mk_upcase_string_with_length(char *s,int n)
{
  char *result = mk_copy_string(s);
  int i;
  for ( i = 0 ; i < n ; i++ )
    if ( result[i] >= 'a' && result[i] <= 'z' )
      result[i] = (char)(result[i] - 'a' + 'A');
  return(result);
}

bool caseless_eq_string(char *s1,char *s2)
{
  int i;
  bool result = TRUE;
  for ( i = 0 ; s1[i] != '\0' && s2[i] != '\0' && result ; i++ )
  {
    char c1 = s1[i];
    char c2 = s2[i];
    if ( c1 >= 'a' && c1 <= 'z' ) c1 = c1 + 'A' - 'a';
    if ( c2 >= 'a' && c2 <= 'z' ) c2 = c2 + 'A' - 'a';
    result = c1 == c2;
  }
  if ( result && s1[i] != s2[i] )
    result = FALSE;
  return result;
}

bool caseless_eq_string_with_length(char *s1,char *s2,int n)
{
  char *d1 = mk_downcase_string(s1);
  char *d2 = mk_downcase_string(s2);
  bool result = eq_string_with_length(d1,d2,n);
  free_string(d1);
  free_string(d2);
  return(result);
}

int count_occurences_of_char_in_string(char *s,char c){
  char *p;
  int count;
  for(p=strchr(s,c),count=0;p;p=strchr(p+1,c),count++);
  return count;
}

int find_char_in_string(char *s,char c){
  char *p = strchr(s,c);
  return (p)? (p-s):-1;
}

void free_string(char *s)
{
  if ( s == NULL )
    my_error("am_free_string NULL s");
  am_free((char *) s,sizeof(char) * (1 + strlen(s)));
}

void am_free_string(char *s)
{
  free_string(s);
}

char *mk_concat_strings(char *s1, char *s2)
{
  char *s = AM_MALLOC_ARRAY(char,1+strlen(s1)+strlen(s2));
  sprintf(s,"%s%s",s1,s2);
  return s;
}

char *mk_amstr_substring(char *s,int start,int end)
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

bool string_pattern_matches(char *pattern, char *string)
{
  //p_pos - index to pattern and s_pos is index to string
  int p_pos = 0, s_pos = 0; 
  bool prev_ast = 0; //true if the previous char in pattern is *

  for(p_pos = 0; p_pos < (int) strlen(pattern); p_pos++)
  {
    if(pattern[p_pos] == '*') 
    {
      prev_ast = 1;
    }
    else
    {
      if(prev_ast)
      {
        /*if the last char was a *, try to match the whole string*/
        while(string[s_pos] != '\0' && pattern[p_pos] != string[s_pos])
        {
          s_pos++;
        }
      }
      if(pattern[p_pos] == string[s_pos])
        s_pos++;
      else
        return 0;
      prev_ast = 0;
    }
  }
  /*last char in pattern wasnot * and string not reached to the end,
  backtrack the string and pattern to make sure they match*/
  if(!prev_ast && s_pos != (int)strlen(string))
  {
    s_pos = (int) strlen(string) - 1;
    p_pos = (int) strlen(pattern) - 1;
    while(pattern[p_pos] != '*' && pattern[p_pos] == string[s_pos])
    {
//printf("%1s matched with %1s\n", &pattern[p_pos], &string[s_pos]);
      s_pos--;
      p_pos--;
    }
    if(pattern[p_pos] != '*')
      return 0;
  }

  return 1;
}

bool old_string_pattern_matches(char *pattern, char *string)
{
  //p_pos - index to pattern and s_pos is index to string
  int p_pos = 0, s_pos = 0; 
  bool prev_ast = 0; //true if the previous char in pattern is *
  char *term;

  for(p_pos = 0; p_pos < (int) strlen(pattern); p_pos++)
  {
    if(pattern[p_pos] == '*') 
    {
      prev_ast = 1;
    }
    else
    {
    /*copy a set of char without * unto term,
      guarantee to have at least 1 char*/
      int p_end = p_pos;
      while(p_end < (int)strlen(pattern) &&
        pattern[p_end] != '*')
      {
        p_end++;
      }
      term = mk_amstr_substring(pattern,p_pos,p_end);
      p_pos = p_end - 1;

      /*try to match the term with string*/
      if(prev_ast)
      {
        /*case when pattern ends in a letter instead of a * */
        if(pattern[p_end] == '\0')
        {
        /*if the term we are matching is bigger than the length in string
          we are looking at*/
          if((int)strlen(term) > (int)strlen(string) - s_pos)
          {
            free_string(term);
            return 0;
          }
          if(eq_string_with_length(term, &string[(int)strlen(string)-(int)strlen(term)], 
            (int)strlen(term)))
          {
            free_string(term);
            return 1;
          }
          else
          {
            free_string(term);
            return 0;
          }
        }

        /*if the last char was a *, try to match the whole string*/
        while(s_pos + (int)strlen(term) <= (int)strlen(string) &&
          !eq_string_with_length(term, &string[s_pos], (int)strlen(term)))
        {
          s_pos++;
        }
      }
      if(!eq_string_with_length(term, &string[s_pos], (int)strlen(term)))
      {
        free_string(term);
        return 0;
      }
      s_pos += (int)strlen(term);
      free_string(term);
      prev_ast = 0;
    }
  }
  /*last char in pattern wasnot * and string not reached to the end,
  return 0*/
  if(!prev_ast && s_pos != (int)strlen(string))
  {
    return 0;
  }

  return 1;
}
