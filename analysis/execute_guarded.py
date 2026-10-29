"""Execute one reviewed local deterministic step with writes confined to its workspace."""
import argparse, os, runpy, sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--workspace',required=True);p.add_argument('script');a=p.parse_args()
root=Path(a.workspace).resolve();script=Path(a.script).resolve()
assert script.is_relative_to(root/'scripts') and script.is_file()
def confined(value):
 if isinstance(value,int):return
 target=Path(os.fsdecode(value)).resolve()
 if not target.is_relative_to(root):raise RuntimeError('BLOCKED_WRITE_OUTSIDE_REPLAY_WORKSPACE: '+str(target))
def audit(event,args):
 if event=='open':
  file,mode,flags=args
  if (isinstance(mode,str) and any(x in mode for x in 'wax+')) or (isinstance(flags,int) and flags & (os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC|os.O_APPEND)):confined(file)
 elif event in ['os.mkdir','os.remove','os.rmdir']:confined(args[0])
 elif event in ['os.rename','os.replace']:confined(args[0]);confined(args[1])
 elif event in ['socket.connect','socket.connect_ex','socket.getaddrinfo','subprocess.Popen','os.system']:
  raise RuntimeError('BLOCKED_NETWORK_OR_PROCESS_OPERATION: '+event)
sys.dont_write_bytecode=True
sys.addaudithook(audit)
runpy.run_path(str(script),run_name='__main__')
