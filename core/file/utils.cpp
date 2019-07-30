/*
 * Copyright (c) 2008-2018 the MRtrix3 contributors.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, you can obtain one at http://mozilla.org/MPL/2.0/
 *
 * MRtrix3 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * For more details, see http://www.mrtrix.org/
 */


#ifdef MRTRIX_MACOSX
# include <sys/param.h>
# include <sys/mount.h>
#elif !defined(MRTRIX_WINDOWS)
# include <sys/vfs.h>
#endif

#include "file/utils.h"


namespace MR
{
  namespace File
  {



    bool use_delayed_writeback (const std::string& path)
    {
#ifdef MRTRIX_WINDOWS
      const unsigned int length = 255;
      char root_path[length];
      if (GetVolumePathName (path.c_str(), root_path, length)) { // Returns non-zero on success

        const unsigned int code = GetDriveType (root_path);
        switch (code) {
          case 0: // DRIVE_UNKNOWN
            DEBUG ("cannot get filesystem information on file \"" + path + "\": " + strerror (errno));
            return true;
          case 1: // DRIVE_NO_ROOT_DIR:
            DEBUG ("erroneous root path derived for file \"" + path + "\": " + strerror (errno));
            return true;
          case 2: // DRIVE_REMOVABLE
            DEBUG ("Drive for file \"" + path + "\" detected as removable; no delayed write-back");
            return false;
          case 3: // DRIVE_FIXED
            DEBUG ("Drive for file \"" + path + "\" detected as fixed; no delayed write-back");
            return false;
          case 4: // DRIVE_REMOTE
            DEBUG ("Drive for file \"" + path + "\" detected as network - using delayed write-back");
            return true;
          case 5: // DRIVE_CDROM
            DEBUG ("Drive for file \"" + path + "\" detected as CD-ROM - using delayed write-back");
            return true;
          case 6: // DRIVE_RAMDISK
            DEBUG ("Drive for file \"" + path + "\" detected as RAM - no delayed write-back");
            return false;
        }

      } else {
        DEBUG ("unable to query root drive path for file \"" + path + "\"; using delayed write-back");
        return true;
      }
#else
      struct statfs fsbuf;
      if (statfs (path.c_str(), &fsbuf)) {
        DEBUG ("cannot get filesystem information on file \"" + path + "\": " + strerror (errno));
        DEBUG ("  defaulting to delayed write-back");
        return true;
      }

      if (fsbuf.f_type == 0xff534d42 /* CIFS */|| fsbuf.f_type == 0x6969 /* NFS */ ||
          fsbuf.f_type == 0x65735546 /* FUSE */ || fsbuf.f_type == 0x517b /* SMB */ ||
          fsbuf.f_type == 0x47504653 /* GPFS */ || fsbuf.f_type == 0xbd00bd0 /* LUSTRE */

#ifdef MRTRIX_MACOSX
          || fsbuf.f_type == 0x0017 /* OSXFUSE */
#endif
      ) {
        DEBUG ("\"" + path + "\" appears to reside on a networked filesystem - using delayed write-back");
        return true;
      }
      DEBUG("\"" + path + "\" does not require delayed write-back");
      return false;
#endif
    }



  }
}
