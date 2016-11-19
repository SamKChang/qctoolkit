import qctoolkit as qtk
import paramiko
import pexpect
import os
import shutil
import subprocess as sp

def remoteRun(cmd, status, session):
  qtk.report('submit', status)
  qtk.report('submit-remote-command', cmd)
  ssh_stdin, ssh_stdout, ssh_stderr = \
    session.exec_command(cmd)
  sshout = trimMsg(ssh_stdout)
  ssherr = trimMsg(ssh_stderr)
  if len(sshout) > 0:
    qtk.report('submit-remote-output', sshout)
  if len(ssherr) > 0:
	  qtk.report('submit-remote-error', ssherr)

def submit(inp_list, root, **remote_settings):
  necessary_list = [
    'ip',
    'submission_script',
  ]
  default_dict = {
    'username': None,
    'password': None,
    'flags': None,
    'timeout': 40,
  }
    
  for k, v in default_dict.iteritems():
    exec "%s = %s" % (k, v)
  if len(inp_list) * 5 > 40:
    timeout = len(inp_list) * 5
  if 'password' in remote_settings:
    password = remote_settings['password']
  if 'username' in remote_settings:
    username = remote_settings['username']
  if 'remote_path' not in remote_settings:
    remote_path = './%s' % root
  else:
    remote_path = remote_settings['remote_path']
  if 'timeout' in remote_settings:
    timeout = remote_settings['timeout']
  if 'prefix' in remote_settings:
    prefix = remote_settings['prefix']
  else:
    prefix = ''
  if 'flags' in remote_settings:
    flags = remote_settings['flags']
  if 'threads' not in remote_settings:
    threads = inp_list[0].setting['threads']
  else:
    threads = remote_settings['threads']
    if threads != inp_list[0].setting['threads']:
      qtk.report('submit', 'reset job threads to %d' % threads)
      for inp in inp_list:
        inp.setting['threads'] = threads
  if 'qthreads' not in remote_settings:
    qthreads = threads
  else:
    qthreads = remote_settings['qthreads']
  for s in necessary_list:
    if s not in remote_settings:
      qtk.exit('cluster setting:%s not defined' % s)
    else:
      exec "%s = '%s'" % (s, remote_settings[s])

  if type(inp_list) is not list:
    inp_list = [inp_list]
  program = inp_list[0].setting['program']

  if os.path.exists(root):
    if 'overwrite' in remote_settings \
    and remote_settings['overwrite']:
      qtk.warning("root directory %s exist, overwrite..." % root)
      shutil.rmtree(root)
      cwd = os.getcwd()
      os.makedirs(root)
      os.chdir(root)
      for inp in inp_list:
        inp.write(inp.molecule.name)
      os.chdir(cwd)
    else:
      qtk.warning("root directory %s exist, uploading existing folder"\
                  % root)
  else:
    cwd = os.getcwd()
    os.makedirs(root)
    os.chdir(root)
    for inp in inp_list:
      inp.write(inp.molecule.name)
    os.chdir(cwd)

  if 'compress' not in remote_settings:
    remote_settings['compress'] = False
    if len(inp_list) > 5:
      remote_settings['compress'] = True

  if remote_settings['compress']:
    qtk.report("submit", "compressing input files")
    cmd = 'tar -zcf %s %s' % (root + '.tar.gz', root)
    run = sp.Popen(cmd, shell=True, stdin=sp.PIPE)
    run.stdin.flush()
    run.communicate()
    run.wait()
    rootToSend = root + '.tar.gz'
    remote_dest = remote_path + '.tar.gz'
    qtk.report("submit", "compression completed")
  else:
    rootToSend = root
    remote_dest = remote_path

  paramiko_kwargs = {}
  if username:
    paramiko_kwargs['username'] = username
  if password:
    paramiko_kwargs['password'] = password
  ssh = paramiko.SSHClient()
  ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
  ssh.load_system_host_keys()
  ssh.connect(ip, **paramiko_kwargs)
  ssh_stdin, ssh_stdout, ssh_stderr = \
    ssh.exec_command('ls %s' % remote_path)
  sshout = ssh_stdout.read()
  if len(sshout) > 0:
    if 'overwrite' in remote_settings \
    and remote_settings['overwrite']:
      status = 'remote path %s exists, overwrite...' % remote_path
      cmd = 'rm -r %s' % remote_path
      remoteRun(cmd, status, ssh)
    else:
      qtk.exit('remote path %s exists' % remote_path)

  ssh_newkey = 'Are you sure you want to continue connecting'
  patterns = [ssh_newkey, '[Pp]assword:', pexpect.EOF]
  if username:
    userStr = username + '@'
  else:
    userStr = ''

  cmd = 'scp -qr %s %s%s:%s' % (rootToSend, userStr, ip, remote_dest)
  qtk.report('submit', 'scp input files...')
  qtk.report('submit-remote_command', cmd)

  p = pexpect.spawn(cmd)
  i = p.expect(patterns, timeout=timeout)
  if i == 0:
    qtk.report('submit', 'adding %s to known_hosts' % ip)
    p.sendline('yes')
    i = p.expect(patterns, timeout=timeout)
  if i == 1:
    p.sendline(password)
    i = p.expect(patterns, timeout=timeout)
  if i == 2:
    if not p.before:
      qtk.report('submit', 'scp completed')
    else:
      qtk.warning('scp message: %s' % p.before)

  if remote_settings['compress']:
    status = "decompress remote input files"
    cmd = 'tar xf %s' % rootToSend
    remoteRun(cmd, status, ssh)
    status = "remove remote tar file"
    cmd = 'rm %s' % rootToSend
    remoteRun(cmd, status, ssh)
    qtk.report('submit', 'done')

  exe = qtk.setting.program_dict[program]
  cmd = "%s \"%s\" %s %d %d '%s' %s" % (submission_script, exe,
    remote_path, threads, qthreads, flags, prefix)
  status = 'submitting jobs...'
  remoteRun(cmd, status, ssh)
  ssh.exec_command("echo %s > %s/cmd.log" % (cmd, remote_path))
  qtk.report('submit', 'submission completed')

  ssh.close()

  if 'debug' in remote_settings and remote_settings['debug']:
    pass
  else:
    qtk.report('submit', 'cleanup local files')
    shutil.rmtree(root)
    if os.path.exists(root + '.tar.gz'):
      os.remove(root + '.tar.gz')

def trimMsg(msg_in):
  out = msg_in.readlines()
  if len(out) > 10:
    out = out[:10]
  return ''.join(filter(None, out))
