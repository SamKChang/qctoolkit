import qctoolkit as qtk
import sqlalchemy as q
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import or_
from datetime import datetime as dt
import datetime
import os

Base = declarative_base()

class Entry(Base):
  __tablename__ = 'entries'

  id = q.Column(q.Integer, primary_key=True)
  date = q.Column(q.DateTime, nullable=False)
  content = q.Column(q.Text)
  comment = q.Column(q.Text)
  data = q.Column(q.Float)

  def __repr__(self):
    if type(self.data) is float:
      return "%s %s %f %s" % (
        self.date.strftime("%Y%m%d-%H:%M:%S"), 
        self.content, 
        self.data,
        self.comment
      )
    else:
      return "%s %s Null %s" % (
        self.date.strftime("%Y%m%d-%H:%M:%S"), 
        self.content, 
        self.comment
      )

class Logger(object):
  def __init__(self, path=':memory:', db_str = None, **kwargs):

    if not db_str:
      db_str = 'sqlite:///' + path

    self.engine = create_engine(db_str, **kwargs)
    self.name = db_str

    if os.path.exists(path):
      qtk.progress('DB', 'loading existing database: %s' % path)

    Base.metadata.create_all(self.engine)
    self.session = self.get_session(new=True)

  def __repr__(self):
    entries = self.list()
    if len(entries) < 5:
      return '%s\n%s' % (self.name, str(entries))
    else:
      return '%s\n%s' % (self.name, str(entries[:5]))

  def push(self, content=None, data=None, comment=None, date=None):
    if date is None:
      date = dt.now()
    entry = Entry(
      date=date, content=content, comment=comment, data=data
    )

    qtk.progress("DB", "adding entry")
    self.session.add(entry)
    try:
      qtk.progress("DB", "attempt to commit...")
      self.session.commit()
      qtk.progress("DB", "done")
    except Exception as err:
      qtk.warning('can not commit, error: %s' % err)

  def list(self, 
    content=None, data=None, comment=None, date=None,
    match = True, epsilon=0.0, dt = datetime.timedelta(0),
    order=False, limit=False, get_list=True):

    if content:
      content_flag = r'%' + content + r'%'
    else:
      content_flag = r'%%'
    if comment:
      comment_flag = r'%' + comment + r'%'
    else:
      comment_flag = r'%%'

    out = self.all(get_list = False)

    if date is not None:
      out = out.filter(
        Entry.date <= date + dt,
        Entry.date >= date - dt,
      )

    if data is not None:
      out = out.filter(
        Entry.data <= data + epsilon,
        Entry.data >= data - epsilon,
      )

    if content is not None:
      if match:
        out = out.filter(Entry.content == content)
      else:
        out = out.filter(Entry.content.like(content_flag))

    if comment is not None:
      if match:
        out = out.filter(Entry.comment == comment)
      else:
        out = out.filter(Entry.comment.like(comment_flag))

    if order == 'ascend':
      out = out.order_by(Entry.data)
    elif order == 'descend':
      out = out.order_by(Entry.data.desc())

    if limit:
      out = out.limit(limit)

    if get_list:
      return out.all()
    else:
      return out

  def all(self, get_list = True):
    if get_list:
      return self.session.query(Entry).all()
    else:
      return self.session.query(Entry)

  def commit(self):
    self.session.commit()

  def get_session(self, new=False):
    if not new:
      return self.session
    else:
      session = sessionmaker(bind=self.engine)
      return session()
