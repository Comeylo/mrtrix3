/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by J-Donald Tournier, 27/06/08.

    This file is part of MRtrix.

    MRtrix is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MRtrix is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "args.h"
#include "app.h"
#include "debug.h"

#define HELP_WIDTH  80

#define HELP_PURPOSE_INDENT 0, 4
#define HELP_ARG_INDENT 8, 20
#define HELP_OPTION_INDENT 2, 20

namespace MR
{
  namespace App
  {

    namespace
    {

      inline size_t size (const std::string& text)
      {
        return text.size() - 2*std::count (text.begin(), text.end(), 0x08U);
      }

      inline void resize (std::string& text, size_t new_size, char fill)
      {
        text.resize (text.size() + new_size - size(text), fill);
      }



      std::string paragraph (
          const std::string& header,
          const std::string& text,
          size_t header_indent,
          size_t indent)
      {
        std::string out, line = std::string (header_indent, ' ') + header + " ";
        if (size (line) < indent)
          resize (line, indent, ' ');

        std::vector<std::string> paragraphs = split (text, "\n");

        for (size_t n = 0; n < paragraphs.size(); ++n) {
          size_t i = 0;
          std::vector<std::string> words = split (paragraphs[n]);
          while (i < words.size()) {
            do {
              line += " " + words[i++];
              if (i >= words.size())
                break;
            }
            while (size (line) + 1 + size (words[i]) < HELP_WIDTH);
            out += line + "\n";
            line = std::string (indent, ' ');
          }
        }
        return out;
      }



      std::string bold (const std::string& text)
      {
        std::string retval (3*text.size(), '\0');
        for (size_t n = 0; n < text.size(); ++n) {
          retval[3*n] = retval[3*n+2] = text[n];
          retval[3*n+1] = 0x08U;
        }
        return retval;
      }


      std::string underline (const std::string& text)
      {
        std::string retval (3*text.size(), '\0');
        for (size_t n = 0; n < text.size(); ++n) {
          retval[3*n] = '_';
          retval[3*n+1] = 0x08U;
          retval[3*n+2] = text[n];
        }
        return retval;
      }

    }





    const char* argtype_description (ArgType type)
    {
      switch (type) {
        case Integer:
          return ("integer");
        case Float:
          return ("float");
        case Text:
          return ("string");
        case ArgFile:
          return ("file");
        case ImageIn:
          return ("image in");
        case ImageOut:
          return ("image out");
        case Choice:
          return ("choice");
        case IntSeq:
          return ("int seq");
        case FloatSeq:
          return ("float seq");
        default:
          return ("undefined");
      }
    }



    std::string help_head (int format) 
    {
      if (!format) 
        return std::string (NAME) + ": part of the MRtrix package\n\n";

      std::string version = MR::printf ("MRtrix %zu.%zu.%zu", VERSION[0], VERSION[1], VERSION[2]);
      std::string date (__DATE__);

      std::string topline = version + std::string (40-size(version)-size(App::NAME)/2, ' ') + bold (App::NAME);
      topline += std::string (80-size(topline)-size(date), ' ') + date;
      
      return topline + "\n\n     " + ( format ? bold (NAME) : std::string (NAME) ) + ": part of the MRtrix package\n\n";
    }




    std::string help_tail (int format) 
    {
      std::string retval;
      if (!format) 
        return retval;

      return bold ("AUTHOR") + "\n" 
        + paragraph ("", AUTHOR, HELP_PURPOSE_INDENT) + "\n"
        + bold ("COPYRIGHT") + "\n" 
        + paragraph ("", COPYRIGHT, HELP_PURPOSE_INDENT) + "\n";
    }







    std::string Description::syntax (int format) const
    {
      std::string s;
      if (format) 
        s += bold ("DESCRIPTION") + "\n\n";
      for (size_t i = 0; i < size(); ++i) 
        s += paragraph ("", (*this)[i], HELP_PURPOSE_INDENT) + "\n";
      return s;
    }




    std::string help_syntax (int format)
    {
      std::string s = "SYNOPSIS";
      if (format)
        s = bold (s) + "\n\n     ";
      else 
        s += ": ";
      s += ( format ? underline (NAME) : NAME ) + " [ options ]";

      for (size_t i = 0; i < ARGUMENTS.size(); ++i) {

        if (ARGUMENTS[i].flags & Optional)
          s += "[";
        s += std::string(" ") + ARGUMENTS[i].id;

        if (ARGUMENTS[i].flags & AllowMultiple) {
          if (! (ARGUMENTS[i].flags & Optional))
            s += std::string(" [ ") + ARGUMENTS[i].id;
          s += " ...";
        }
        if (ARGUMENTS[i].flags & (Optional | AllowMultiple))
          s += " ]";
      }
        return s + "\n\n";
      }




    std::string Argument::syntax (int format) const
    {
      std::string retval = paragraph (( format ? underline (id) : id ), desc, HELP_ARG_INDENT);
      if (format) 
        retval += "\n";
      return retval;
    }





    std::string ArgumentList::syntax (int format) const
    {
      std::string s;
      for (size_t i = 0; i < size(); ++i)
        s += (*this)[i].syntax (format);
      return s + "\n";
    }





    std::string Option::syntax (int format) const
    {
      std::string opt ("-");
      opt += id;

      if (format)
        opt = underline (opt);

      for (size_t i = 0; i < size(); ++i)
        opt += std::string (" ") + (*this)[i].id;

      if (format) 
        opt = "  " + opt + "\n" + paragraph ("", desc, HELP_PURPOSE_INDENT);
      else
        opt = paragraph (opt, desc, HELP_OPTION_INDENT);
      if (format) 
        opt += "\n";
      return opt;
    }




    std::string OptionGroup::syntax (int format) const
    {
      std::string s ( format ? bold (name) + "\n\n" : std::string (name) + ":\n" );

      for (size_t i = 0; i < size(); ++i)
        s += (*this) [i].syntax (format);
      if (!format)
        s += "\n";
      return s;
    }




    std::string OptionList::syntax (int format) const
    {
      std::string s;
      for (size_t i = 0; i < size(); ++i)
        s += (*this) [i].syntax (format);
      return s;
    }








    std::string Argument::usage () const
    {
      std::ostringstream stream;
      stream << "ARGUMENT " << id << " "
        << (flags & Optional ? '1' : '0') << " "
        << (flags & AllowMultiple ? '1' : '0') << " ";

      switch (type) {
        case Integer:
          stream << "INT " << defaults.i.min << " " << defaults.i.max << " " << defaults.i.def;
          break;
        case Float:
          stream << "FLOAT " << defaults.f.min << " " << defaults.f.max << " " << defaults.f.def;
          break;
        case Text:
          stream << "TEXT";
          if (defaults.text)
            stream << " " << defaults.text;
          break;
        case ArgFile:
          stream << "FILE";
          break;
        case Choice:
          stream << "CHOICE";
          for (const char* const* p = defaults.choices.list; *p; ++p)
            stream << " " << *p;
          stream << " " << defaults.choices.def;
          break;
        case ImageIn:
          stream << "IMAGEIN";
          break;
        case ImageOut:
          stream << "IMAGEOUT";
          break;
        case IntSeq:
          stream << "ISEQ";
          break;
        case FloatSeq:
          stream << "FSEQ";
          break;
        default:
          assert (0);
      }
      stream << "\n";
      if (desc.size())
        stream << desc << "\n";

      return stream.str();
    }




    std::string Option::usage () const
    {
      std::ostringstream stream;
      stream << "OPTION " << id << " "
        << (flags & Optional ? '1' : '0') << " "
        << (flags & AllowMultiple ? '1' : '0') << "\n";

      if (desc.size())
        stream << desc << "\n";

      for (size_t i = 0; i < size(); ++i)
        stream << (*this)[i].usage ();
      
      return stream.str();
    }


  }
}


