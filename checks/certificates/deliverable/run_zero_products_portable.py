#!/usr/bin/env python3
"""Lancer le verificateur original des zeros avec mesure memoire sous Windows.

Le code mathematique de check_zero_products.py reste inchange et son empreinte
est verifiee dans SOURCE-MANIFEST.json avant execution. Sous Windows, le seul
service resource utilise (ru_maxrss en Kio) est fourni par PeakWorkingSetSize
du processus courant. La limite de 900000 Kio et les controles de temps de
l'original restent actifs. Aucune valeur de memoire n'est simulee.
API Windows : https://learn.microsoft.com/en-us/windows/win32/api/psapi/ns-psapi-process_memory_counters
"""
import argparse, ctypes, hashlib, json, os, runpy, sys, types
from pathlib import Path

def install_windows_resource():
    from ctypes import wintypes
    class ProcessMemoryCounters(ctypes.Structure):
        _fields_=[('cb',wintypes.DWORD),('PageFaultCount',wintypes.DWORD),
          ('PeakWorkingSetSize',ctypes.c_size_t),('WorkingSetSize',ctypes.c_size_t),
          ('QuotaPeakPagedPoolUsage',ctypes.c_size_t),('QuotaPagedPoolUsage',ctypes.c_size_t),
          ('QuotaPeakNonPagedPoolUsage',ctypes.c_size_t),('QuotaNonPagedPoolUsage',ctypes.c_size_t),
          ('PagefileUsage',ctypes.c_size_t),('PeakPagefileUsage',ctypes.c_size_t)]
    kernel=ctypes.WinDLL('kernel32',use_last_error=True)
    psapi=ctypes.WinDLL('psapi',use_last_error=True)
    kernel.GetCurrentProcess.restype=wintypes.HANDLE
    psapi.GetProcessMemoryInfo.argtypes=[wintypes.HANDLE,ctypes.POINTER(ProcessMemoryCounters),wintypes.DWORD]
    psapi.GetProcessMemoryInfo.restype=wintypes.BOOL
    def getrusage(who):
        if who!=0:raise ValueError('Seul RUSAGE_SELF est pris en charge')
        data=ProcessMemoryCounters();data.cb=ctypes.sizeof(data)
        if not psapi.GetProcessMemoryInfo(kernel.GetCurrentProcess(),ctypes.byref(data),data.cb):
            raise ctypes.WinError(ctypes.get_last_error())
        return types.SimpleNamespace(ru_maxrss=(int(data.PeakWorkingSetSize)+1023)//1024)
    adapter=types.ModuleType('resource');adapter.RUSAGE_SELF=0;adapter.getrusage=getrusage
    sys.modules['resource']=adapter
    print('Windows memory metric: actual PeakWorkingSetSize, rounded up to KiB.',flush=True)

def main():
    p=argparse.ArgumentParser(add_help=False);p.add_argument('--work',required=True,type=Path);p.add_argument('--output',required=True,type=Path)
    a,_=p.parse_known_args();work=a.work.resolve();output=a.output.resolve()
    manifest=json.loads((work/'SOURCE-MANIFEST.json').read_text())
    row=next(r for r in manifest['new_verifiers'] if r['name']=='check_zero_products.py')
    script=work/row['name'];data=script.read_bytes()
    if len(data)!=row['bytes'] or hashlib.sha256(data).hexdigest()!=row['sha256']:
        raise RuntimeError('Le verificateur original des zeros a change')
    if os.name=='nt':install_windows_resource()
    runpy.run_path(str(script),run_name='__main__')
    # Canonical serialization only: Windows path separators are not mathematical labels.
    result=json.loads(output.read_text())
    hashes={key.replace('\\','/'):value for key,value in result['manifest_hashes'].items()}
    if len(hashes)!=len(result['manifest_hashes']):raise RuntimeError('Collision de chemins')
    result['manifest_hashes']=hashes
    for record in result['records']:record['manifest']=record['manifest'].replace('\\','/')
    result['portable_paths']='POSIX separators; mathematical labels and results unchanged'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print('Portable receipt: canonical POSIX path separators.',flush=True)
if __name__=='__main__':main()
